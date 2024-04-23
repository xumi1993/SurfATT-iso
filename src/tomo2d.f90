!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: Jan 2024, Initialize Codes
!
!=====================================================================
module tomo2d

  use para, ap => att_para_global
  use src_rec, sr_ph => src_rec_global_ph, sr_gr => src_rec_global_gr
  use grid, ag_ph => att_grid_global_ph, ag_gr => att_grid_global_gr
  use measadj
  use acqui_2d
  use setup_att_log

  implicit none

  type, public :: att_tomo_2d
    character(len=MAX_NAME_LEN)                            :: module = "TOMO2D"
    character(len=MAX_STRING_LEN)                          :: message
    contains
    procedure :: init => init_tomo_2d, do_forward, initialize_inv, do_inversion, &
                 post_proc_eikokernel, linesearch, optimize, forward_simulate,break_iter,&
                 steepest_descent
  end type

  type(att_acqui_2d), pointer, private                     :: acqui
  integer                                                  :: itype=1

contains

  subroutine init_tomo_2d(this)
    class(att_tomo_2d), intent(inout) :: this
    character(MAX_STRING_LEN) :: errmsg
    integer :: stat, i, j
    integer, dimension(:,:), allocatable :: isrcs

    do itype = 1, 2
      if (.not. ap%data%vel_type(itype)) cycle
      call select_type()
      call acqui%init(itype)
      call acqui%sr%scatter_src_gather()
    enddo
    if (myrank == 0) then
      call EXECUTE_COMMAND_LINE('mkdir -p '//trim(ap%output%output_path),&
                                exitstat=stat, cmdmsg=errmsg)
      if (stat /= 0) then
        call write_log(errmsg, 3, this%module)
        call exit_MPI(myrank, errmsg)
      endif
      if (ap%inversion%sigma_2d == 0.0_dp) then
        call write_log('Please setup sigma_2d for smoothing.',3,this%module)
        stop
      endif
    endif
    call synchronize_all()
  end subroutine init_tomo_2d

  subroutine do_forward(this, ncb, pert_vel, hmarg)
    class(att_tomo_2d), intent(inout) :: this
    integer, intent(in) :: ncb(2)
    real(kind=dp), intent(in) :: pert_vel, hmarg
    type(att_measadj) :: ma
    integer :: i
    character(len=MAX_STRING_LEN) :: fname
    real(kind=dp), dimension(:), allocatable :: global_tt
    real(kind=dp) :: chi_global

    ! loop for sources
    do itype = 1, 2
      if (.not. ap%data%vel_type(itype)) cycle
      call select_type()
      if ((ncb(1) > 0) .and. (ncb(2) > 0)) then
        call acqui%add_pert(ncb(1), ncb(2), pert_vel, hmarg)
        call acqui%write_target_model()
      endif
      call acqui%prepare_fwd()
      call write_log("This is forward simulation for "//ap%data%gr_name(itype),1,this%module)
      call this%forward_simulate(chi_global, .true., .false.)
      write(fname, '(a,"/src_rec_file_forward_",a,".csv")') &
          trim(ap%output%output_path),trim(ap%data%gr_name(itype))
      call acqui%sr%to_csv(fname)
    enddo
    call synchronize_all()
  end subroutine do_forward

  subroutine do_inversion(this)
    class(att_tomo_2d), intent(inout) :: this
    real(kind=dp) :: chi_global
    character(len=MAX_STRING_LEN) :: fname
    logical :: isbreak
    integer :: i, it
    
    do itype = 1, 2
      if (.not. ap%data%vel_type(itype)) cycle
      call select_type()
      call this%initialize_inv()
      do it = 1, ap%inversion%niter 
        acqui%iter = it
        write(this%message, '(a,i0,a)') 'Iteration: ',acqui%iter, 'th'
        call write_log(this%message,1,this%module)
        ! Prepare velocity for forward simulation
        call acqui%prepare_fwd()
        ! Prepare kernel matrix for adjoint simulation
        call acqui%prepare_inv()
        ! forward and adjoint simulation
        call this%forward_simulate(chi_global, .true., .true.)
        call acqui%write_obj_func()
        ! write objective function
        if (myrank==0) then
          if (acqui%iter == 1) acqui%chi0 = chi_global
          acqui%misfits(acqui%iter) = chi_global
          write(this%message, '(a,F0.4," (",F0.4,"%)")') 'Total misfit of '//&
                trim(ap%data%gr_name(itype))//': ',&
                chi_global,100*chi_global/acqui%chi0
          ! write synthetic tt to file
          write(fname, '(a,I3.3,".csv")') trim(ap%output%output_path)&
              //'src_rec_file_forward_'//trim(ap%data%gr_name(itype))//'_',acqui%iter-1
          call write_log(this%message,1,this%module)
          call this%break_iter(isbreak)
        endif ! myrank==0
        call bcast_all(isbreak)
        if (isbreak) exit
        call acqui%sr%to_csv(trim(fname))
        ! optimization
        call this%optimize()
        call acqui%write_iter()
      enddo
      ! write final model
      call acqui%write_model()
      if(myrank == 0) call acqui%h%close()
      call write_log('Inversion is done.',1,this%module)
    enddo
  end subroutine do_inversion

  subroutine initialize_inv(this)
    class(att_tomo_2d), intent(inout) :: this
    if (myrank == 0) then 
      open(IOUT, file=trim(ap%output%output_path)//'/objective_function_'//&
           trim(ap%data%gr_name(itype))//'.csv', status='replace',action='write')
      call acqui%h%open(acqui%model_fname, status='new', action='write')    
    endif
    acqui%updatemax = ap%inversion%step_length
    call acqui%write_iter()
    call synchronize_all()
  end subroutine initialize_inv

  subroutine post_proc_eikokernel(this, pidx, adjtable, timetable)
    class(att_tomo_2d), intent(inout) :: this
    integer, intent(inout) :: pidx
    real(kind=dp), dimension(:,:),allocatable :: Tx, Ty
    real(kind=dp), dimension(:,:,:), allocatable :: adj_s_local
    real(kind=dp),  dimension(:,:), allocatable, intent(in) :: adjtable, timetable
    ! real(kind=dp), dimension(acqui%ag%nx, acqui%ag%ny) :: vtmp
    integer :: i, j, k

    ! vtmp = 1/acqui%ag%svel(pidx, :,:)
    acqui%adj_s_local(pidx, :,:) = acqui%adj_s_local(pidx, :,:)+adjtable / acqui%ag%svel(pidx, :,:)**2
  end subroutine post_proc_eikokernel

  subroutine forward_simulate(this, chi, istotable, isadj)
    class(att_tomo_2d), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp) :: chi_local, chi
    real(kind=dp), dimension(:,:), allocatable :: adj
    real(kind=dp), dimension(:), allocatable :: local_tt
    integer :: i, j
    logical :: istotable, isadj

    local_tt = zeros(acqui%sr%npath)
    if ((acqui%sr%iend-acqui%sr%istart)>=0) then
      chi_local = 0.0_dp
      do i = acqui%sr%istart, acqui%sr%iend
        ! write log
        write(this%message, '(a,F0.4,a,a)') 'period: ',&
              acqui%sr%periods(acqui%sr%isrcs(i,1)),' src_name: ',&
              acqui%sr%stations%staname(acqui%sr%isrcs(i, 2))
        call write_log(this%message,0,this%module)
        ! get receiver gather
        call ma%get_recs(acqui%sr,acqui%sr%isrcs(i,1),acqui%sr%stations%staname(acqui%sr%isrcs(i, 2)))
        ! forward
        call ma%run_forward(acqui%ag%svel(acqui%sr%isrcs(i,1),:,:), acqui%ag%m11(acqui%sr%isrcs(i, 1),:,:),&
                            acqui%ag%m22(acqui%sr%isrcs(i,1),:,:), acqui%ag%m12(acqui%sr%isrcs(i,1),:,:),&
                            acqui%ag%ref_t(acqui%sr%isrcs(i,1),:,:))
        ! to time table
        if (istotable) call ma%to_table(local_tt)
        ! measure adjoint
        if (isadj) then
          call ma%run_adjoint(acqui%ag%m11(acqui%sr%isrcs(i, 1),:,:),acqui%ag%m22(acqui%sr%isrcs(i,1),:,:),&
                              acqui%ag%m12(acqui%sr%isrcs(i,1),:,:),adj)
          ! post proc of eikonal kernel
          call this%post_proc_eikokernel(acqui%sr%isrcs(i,1), adj, ma%timetable)
        endif
        ! sum chi
        chi_local = chi_local + ma%chi
        call ma%distory()
      enddo ! i = acqui%istart, acqui%iend
    endif ! if ((acqui%iend-acqui%istart)>=0)
    call synchronize_all()
    call sum_all(acqui%adj_s_local, acqui%adj_s, acqui%sr%nperiod,acqui%ag%nx,acqui%ag%ny)
    call sum_all(chi_local, chi)
    call bcast_all(chi)
    if (istotable) call sum_all(local_tt, acqui%sr%tt_fwd,acqui%sr%npath)
  end subroutine forward_simulate

  subroutine optimize(this)
    class(att_tomo_2d), intent(inout) :: this
    integer :: istart, iend, i
    real(kind=dp) :: sigma
    real(kind=dp), dimension(:,:), allocatable :: tmp
    real(kind=dp), dimension(:,:,:), allocatable :: ker_s_local
    
    call scatter_all_i(acqui%ag%nperiod,mysize,myrank,istart,iend)
    call write_log('This is optimization...',1,this%module)
    ker_s_local = zeros(acqui%sr%nperiod, acqui%ag%nx, acqui%ag%ny)
    ! acqui%ker_s = 0.0_dp
    ! regularization
    if (iend-istart >= 0) then
      do i = istart, iend
        write(this%message, '(a,F0.4,a)') 'Optimization for period: ',acqui%ag%periods(i),'s' 
        call write_log(this%message,0,this%module)
        ! sigma = km2deg*ap%topo%wavelen_factor*acqui%sr%periods(i)*sum(acqui%svel(i,:,:))/size(acqui%svel(i,:,:))
        sigma = ap%inversion%sigma_2d
        ker_s_local(i,:,:)  = gaussian_smooth_geo_2(acqui%adj_s(i,:,:),acqui%ag%xgrids,acqui%ag%ygrids,sigma)
      enddo
    endif
    call synchronize_all()
    call sum_all(ker_s_local, acqui%ker_s,acqui%sr%nperiod,acqui%ag%nx,acqui%ag%ny)
    if (ap%inversion%optim_method == 0) then
      ! steepest descent
      call this%steepest_descent()
    else
      ! line search
      call this%linesearch()
    endif
  end subroutine optimize

  subroutine steepest_descent(this)
    class(att_tomo_2d), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), allocatable :: gradient
  
    call write_log('Optimizing using steepest descent...',1,this%module)
    if (acqui%misfits(acqui%iter) > acqui%misfits(acqui%iter-1)) then
      write(this%message, '(a,f0.8,a,f0.8)') 'Misfit increased from ',&
            acqui%misfits(acqui%iter-1),' to ',acqui%misfits(acqui%iter)
      call write_log(this%message,1,this%module)
      acqui%updatemax = acqui%updatemax*ap%inversion%maxshrink
      write(this%message, '(a,F0.4)') 'Shrink step length to ',acqui%updatemax
      call write_log(this%message, 1, this%module)
    endif
    if (local_rank == 0) then
      gradient = acqui%ker_s/maxval(abs(acqui%ker_s))
      acqui%svel = acqui%svel*(1+acqui%updatemax * gradient)
    endif
    call synchronize_all()

  end subroutine steepest_descent

  subroutine linesearch(this)
    class(att_tomo_2d), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp), dimension(:,:,:), allocatable :: update_s,svel_new
    real(kind=dp) :: chi0, chi_local, chi_global, chi, qpt
    real(kind=dp), dimension(:,:), allocatable :: adj
    integer :: i, j, sit

    call write_log('This is line search...',1,this%module)
    ! initialize
    chi0 = acqui%misfits(acqui%iter)
    ! do line search
    do sit = 1, ap%inversion%max_sub_niter
      write(this%message, '(a,i0,a,f6.4)') 'Sub-iteration: ',sit, 'th, step length: ',acqui%updatemax
      call write_log(this%message,1,this%module)
      call acqui%prepare_fwd_linesearch()
      call this%forward_simulate(chi_global, .false., .false.)
      call bcast_all(chi_global)
      if (chi_global > chi0) then
        write(this%message, '(a,f6.4,a,f6.4)') 'Misfit increased from ',chi0,' to ',chi_global
        call write_log(this%message,0,this%module)
        acqui%updatemax = acqui%updatemax*ap%inversion%maxshrink
      else
        write(this%message, '(a,f6.4,a,f6.4)') 'Misfit decreased from ',chi0,' to ',chi_global
        call write_log(this%message, 0,this%module)
        exit
      endif
    enddo
    if (local_rank == 0) then
      acqui%svel = acqui%ag%svel
    endif
    call synchronize_all()    
  end subroutine linesearch

  subroutine break_iter(this, isbreak)
    class(att_tomo_2d), intent(inout) :: this
    real(kind=dp) :: misfit_prev, misfit_curr, misfit_diff
    logical, intent(out) :: isbreak

    isbreak = .false.
    if (acqui%iter > iter_store) then
      misfit_prev = sum(acqui%misfits(acqui%iter-iter_store:acqui%iter-1))/iter_store
      misfit_curr = sum(acqui%misfits(acqui%iter-iter_store+1:acqui%iter))/iter_store
      misfit_diff = (misfit_prev-misfit_curr)/misfit_prev
      write(this%message,'(a,i2,a,F10.8)') 'Misfit change of last', &
            iter_store,' iterations: ',misfit_diff
      call write_log(this%message,1,this%module)
      if (abs(misfit_diff) < ap%inversion%min_derr) isbreak = .true.
    endif    
  end subroutine break_iter

  subroutine select_type()
    if (itype == 1) then
      acqui => acqui_ph
      ap%data%igr = 0
    else
      acqui => acqui_gr
      ap%data%igr = 1
    endif
  end subroutine select_type

end module