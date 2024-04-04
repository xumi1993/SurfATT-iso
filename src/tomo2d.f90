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
                 post_proc_eikokernel, linesearch, optimize, forward_simulate,break_iter
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
      call acqui%scatter_src_gather()
    enddo
    if (myrank == 0) then
      call EXECUTE_COMMAND_LINE('mkdir -p '//trim(ap%output%output_path),&
                                exitstat=stat, cmdmsg=errmsg)
      if (stat /= 0) then
        call write_log(errmsg, 3, this%module)
        call exit_MPI(myrank, errmsg)
      endif
      ! call initialize_files()
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
      if ((ncb(1) > 0) .and. (ncb(2) > 0)) call acqui%add_pert(ncb(1), ncb(2), pert_vel, hmarg)
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
        ! write objective function
        if (myrank==0) then
          if (acqui%iter == 1) acqui%chi0 = chi_global
          acqui%misfits(acqui%iter) = chi_global
          write(this%message, '(a,F0.2," (",F0.2,"%)")') 'Total misfit of '//&
                trim(ap%data%gr_name(itype))//': ',&
                chi_global,100*chi_global/acqui%chi0
          ! write synthetic tt to file
          write(fname, '(a,I3.3,".csv")') trim(ap%output%output_path)&
              //'src_rec_file_forward_'//trim(ap%data%gr_name(itype))//'_',acqui%iter-1
          call write_log(this%message,1,this%module)
          call this%break_iter(isbreak)
          if(isbreak) exit
        endif ! myrank==0
        call acqui%sr%to_csv(trim(fname))
        ! optimization
        call this%optimize()
        call acqui%write_iter()
      enddo
      ! write final model
      call acqui%write_model()
      call acqui%h%close()
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
    ! call init_model()
    call synchronize_all()
  end subroutine initialize_inv

  subroutine post_proc_eikokernel(this, pidx, adjtable, timetable)
    class(att_tomo_2d), intent(inout) :: this
    integer, intent(inout) :: pidx
    real(kind=dp), dimension(:,:),allocatable :: Tx, Ty
    ! real(kind=dp), dimension(:,:,:,:), allocatable :: adj_s_local, adj_xi_local, adj_eta_local
    real(kind=dp),  dimension(:,:), allocatable, intent(in) :: adjtable, timetable
    real(kind=dp), dimension(acqui%ag%nx, acqui%ag%ny) :: vtmp, lat_corr
    integer :: i, j, k

    vtmp = 1/acqui%ag%svel(pidx, :,:)
    acqui%adj_s(pidx, :,:) = acqui%adj_s(pidx, :,:)+adjtable * vtmp**3
  end subroutine post_proc_eikokernel

  subroutine forward_simulate(this, chi, istotable, isadj)
    class(att_tomo_2d), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp) :: chi_local, chi
    real(kind=dp), dimension(:,:), allocatable :: adj
    real(kind=dp), dimension(acqui%sr%npath) :: local_tt
    integer :: i, j
    logical :: istotable, isadj

    if ((acqui%iend-acqui%istart)>=0) then
      chi_local = 0.0_dp
      do i = acqui%istart, acqui%iend
        ! write log
        write(this%message, '(a,F0.4,a,a)') 'period: ',&
              acqui%sr%periods(acqui%isrcs(i,1)),' src_name: ',&
              acqui%sr%stations%staname(acqui%isrcs(i, 2))
        call write_log(this%message,0,this%module)
        ! get receiver gather
        call ma%get_recs(acqui%sr,acqui%isrcs(i,1),acqui%sr%stations%staname(acqui%isrcs(i, 2)))
        ! forward
        call ma%run_forward(acqui%ag%svel(acqui%isrcs(i,1),:,:), acqui%ag%m11(acqui%isrcs(i, 1),:,:),&
                            acqui%ag%m22(acqui%isrcs(i,1),:,:), acqui%ag%m12(acqui%isrcs(i,1),:,:),&
                            acqui%ag%ref_t(acqui%isrcs(i,1),:,:))
        ! to time table
        if (istotable) call ma%to_table(local_tt)
        ! measure adjoint
        if (isadj) then
          call ma%run_adjoint(acqui%ag%m11(acqui%isrcs(i, 1),:,:),acqui%ag%m22(acqui%isrcs(i,1),:,:),&
                              acqui%ag%m12(acqui%isrcs(i,1),:,:),adj)
          ! post proc of eikonal kernel
          call this%post_proc_eikokernel(acqui%isrcs(i,1), adj, ma%timetable)
        endif
        ! sum chi
        chi_local = chi_local + ma%chi
        call ma%distory()
      enddo ! i = acqui%istart, acqui%iend
    endif ! if ((acqui%iend-acqui%istart)>=0)
    call synchronize_all()
    call sum_all_dp(chi_local, chi)
    if (istotable) call sum_all_1Darray_dp(local_tt, acqui%sr%tt_fwd,acqui%sr%npath)
  end subroutine forward_simulate

  subroutine optimize(this)
    class(att_tomo_2d), intent(inout) :: this
    integer :: istart, iend, i
    real(kind=dp) :: sigma
    real(kind=dp), dimension(:,:), allocatable :: tmp
    
    call scatter_all_i(acqui%ag%nperiod,mysize,myrank,istart,iend)
    call write_log('This is optimization...',1,this%module)
    acqui%ker_s = 0.0_dp
    ! regularization
    do i = istart, iend
      if (iend-istart >= 0) then
        write(this%message, '(a,F0.4,a)') 'Optimization for period: ',acqui%ag%periods(i),'s' 
        call write_log(this%message,0,this%module)
        sigma = km2deg*ap%topo%wavelen_factor*sum(acqui%svel(i,:,:))/size(acqui%svel(i,:,:))
        acqui%ker_s(i,:,:)  = gaussian_smooth_geo_2(acqui%adj_s(i,:,:),acqui%ag%xgrids,acqui%ag%ygrids,sigma)
        ! acqui%ker_s(i,:,:) = tmp
      endif
    enddo
    call synchronize_all()
    ! line search
    call this%linesearch()
  end subroutine optimize

  subroutine linesearch(this)
    class(att_tomo_2d), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp), dimension(:,:,:), allocatable :: update_s, update_xi, update_eta,&
                                                    svel_new, xi_new, eta_new
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
      call acqui%prepare_fwd_linesearch(xi_new, eta_new)
      call this%forward_simulate(chi_global, .false., .false.)
      call bcast_all_singledp(chi_global)
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
    if (myrank == 0) then
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
      write(this%message,'(a,i2,a,F6.4)') 'Misfit change of last', &
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