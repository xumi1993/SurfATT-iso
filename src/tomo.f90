!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: Oct 2023, Initialize Codes
!
!=====================================================================
module tomo
  use para, ap => att_para_global
  use model, am => att_model_global
  use acqui, aq_ph => att_acqui_global_ph, aq_gr => att_acqui_global_gr
  use decomposer, amd => att_mesh_decomposer_global
  use measadj
  use utils
  use hdf5_interface
  use surfker, only: depthkernel_mpi
  use setup_att_log

  implicit none

  integer, private :: iter, itype, OID
  real(kind=dp), private :: updatemax
  character(len=MAX_STRING_LEN) :: model_fname
  type(hdf5_file) :: h

  type, public ::  att_tomo
    integer, dimension(:,:), allocatable                   :: isrcs
    integer                                                :: nsrc, istart, iend,iter
    character(len=MAX_STRING_LEN)                          :: message,module='TOMO'
    real(kind=dp), dimension(:), pointer                   :: misfits
    logical                                                :: is_fwd
  contains
    procedure :: init => initialize_tomo
    procedure :: do_forward, do_inversion
    procedure, private :: eikokernel,break_iter,initialize_inv, line_search, reset_opt, steepest_descent
  end type

  integer :: win_misfits
  type(att_acqui), pointer :: aq

contains

  subroutine initialize_tomo(this, is_fwd)
    class(att_tomo), intent(inout) :: this
    logical, optional, intent(in) :: is_fwd
    integer :: stat
    character(len=MAX_STRING_LEN) :: errmsg

    if (present(is_fwd)) then
      this%is_fwd = is_fwd
    else
      this%is_fwd = .false.
    endif

    if (myrank == 0) then
      call EXECUTE_COMMAND_LINE('mkdir -p '//trim(ap%output%output_path),&
                                exitstat=stat, cmdmsg=errmsg)
      if (stat /= 0) then
        call write_log(errmsg, 3, this%module)
        call exit_MPI(myrank, errmsg)
      endif
      model_fname = trim(ap%output%output_path)//"/"//modfile
      ! call initialize_files()
    endif
    do itype = 1, 2
      if(.not. ap%data%vel_type(itype)) cycle
      call select_type()
      call aq%init(itype, this%is_fwd)
    enddo
    call am%prepare_model_opt()
    call synchronize_all()
  end subroutine initialize_tomo

  subroutine initialize_inv(this)
    class(att_tomo), intent(inout) :: this
    
    updatemax = ap%inversion%step_length
    call prepare_shm_array_dp_1d(this%misfits, ap%inversion%niter, win_misfits)
    if (myrank == 0) then 
      open(OID, file=trim(ap%output%output_path)//'/objective_function.csv',&
          status='replace',action='write')
      if (ap%output%verbose_level > 0) then
        call h%open(model_fname, status='new', action='w')
        do itype = 1, 2
          if (.not. ap%data%vel_type(itype)) cycle
          call select_type()
          call h%add('/stlo_'//trim(ap%data%gr_name(itype)), aq%sr%stations%stlo)
          call h%add('/stla_'//trim(ap%data%gr_name(itype)), aq%sr%stations%stla)
        enddo
        call h%add('/lon', am%xgrids)
        call h%add('/lat', am%ygrids)
        call h%add('/dep', am%zgrids)
        call h%add('/vs_000', transpose_3(am%vs3d))
      endif
    endif
    call synchronize_all()
  end subroutine initialize_inv

  subroutine do_forward(this, max_noise)
    class(att_tomo), intent(inout) :: this
    real(kind=dp), optional, intent(in) :: max_noise
    type(att_measadj) :: ma
    integer :: i
    character(len=MAX_STRING_LEN) :: fname
    real(kind=dp) :: chi_global
    logical :: add_noise

    if (present(max_noise)) then
      add_noise = .true.
    else
      add_noise = .false.
    endif
    call this%reset_opt()
    ! loop for sources
    do itype = 1, 2
      if (.not. ap%data%vel_type(itype)) cycle
      call select_type()
      call aq%forward_simulate(chi_global, .true., .false.)
      if (myrank == 0 .and. add_noise) then
        call aq%sr%add_random_noise(max_noise)
      endif
      write(fname, '(a,"/src_rec_file_forward_",a,".csv")') &
          trim(ap%output%output_path),trim(ap%data%gr_name(itype))
      call aq%sr%to_csv(fname)
    enddo
    call write_log('Forward simulation is done.',1,this%module)
  end subroutine do_forward

  subroutine do_inversion(this)
    class(att_tomo), intent(inout) :: this
    integer :: i, iz, ip
    logical :: status_ok
    real(kind=dp), dimension(:,:,:), allocatable :: update
    character(len=MAX_STRING_LEN) :: fname
    real(kind=dp), dimension(:,:), allocatable :: adj

    call this%initialize_inv()
    do iter = 1, ap%inversion%niter
      call this%reset_opt()
      ! starting simulation of travel time for each src-period pair
      write(this%message, '(a,i0,a)')'Start ',iter-1,'th iteration...'
      call write_log(this%message,1,this%module)
      do itype = 1, 2
        if (.not. ap%data%vel_type(itype)) cycle
        call select_type()
        call aq%prepare_inv()
        ! compute eikonal kernel
        call this%eikokernel()
        ! compute depth kernel
        call aq%depthkernel()
        ! combine tt kernel with depth kernel
        call aq%combine_kernels()
        ! optimzation
        call aq%post_proc()
      enddo
      if (myrank == 0) then
        ! write objective function
        write(OID, '(i3,",",f10.2,",",f6.4)') &
              iter-1,this%misfits(iter),updatemax
        call flush(OID)
      endif
      call this%break_iter(status_ok)
      if (status_ok) exit
      if (ap%inversion%optim_method == 0) then
        call this%steepest_descent()
      else
        call this%line_search()
      endif
      call synchronize_all()
    enddo
    ! write final model
    call am%write('final_model')
    if (myrank == 0 .and. ap%output%verbose_level > 0) call h%close(finalize=.true.)
    close(OID)
    call write_log('Inversion is done.',1,this%module)
  end subroutine do_inversion

  subroutine eikokernel(this)
    ! FILEPATH: /Users/xumijian/Codes/SurfATT/src/tomo.f90
    !> Compute the eikonal kernel for seismic tomography.
    !! Inputs:
    !!   this: an object of type att_tomo
    class(att_tomo), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp) :: chi_global 
    real(kind=dp), dimension(:,:), allocatable :: adj, kden
    character(len=MAX_STRING_LEN) :: fname
    integer :: iz,i

    ! Forward simulation for surface wave velocity
    call aq%forward_simulate(chi_global, .true., .true.)
    if (myrank == 0) then
      if (iter == 1) aq%chi0 = chi_global
      this%misfits(iter) =  this%misfits(iter)+chi_global
      write(this%message, '(a,F0.2," (",F0.2,"%)")') 'Total misfit of '//&
            trim(ap%data%gr_name(itype))//': ',&
            chi_global,100*chi_global/aq%chi0
      call write_log(this%message,1,this%module)
    endif
    ! write synthetic tt to file
    if (ap%output%verbose_level > 1) then
      if (myrank == 0) write(fname, '(a,I3.3,".csv")') trim(ap%output%output_path)&
          //'src_rec_file_forward_'//trim(ap%data%gr_name(itype))//'_',iter-1
      call aq%sr%to_csv(trim(fname))
    endif
    call synchronize_all()
    
  end subroutine eikokernel

  subroutine steepest_descent(this)
    class(att_tomo), intent(inout) :: this
    real(kind=dp) :: chi, chi_local
    real(kind=dp) :: max_gk
    character(len=MAX_STRING_LEN) :: secname

    call write_log('Optimizing using steepest descent...',1,this%module)
    if (myrank == 0) then
      call write_gradient()
      max_gk = maxval(abs(gradient_s))
      if (iter>1 .and. this%misfits(iter) > this%misfits(iter-1)) then
        write(this%message, '(a,F0.4,a,F0.4)') 'Misfit increase from ',this%misfits(iter-1),' to ',this%misfits(iter)
        call write_log(this%message, 1, this%module)
        updatemax = updatemax * ap%inversion%maxshrink
        write(this%message, '(a,F0.4)') 'Shrink step length to ',updatemax
        call write_log(this%message, 1, this%module)
      endif
      gradient_s = -updatemax * gradient_s / max_gk
      am%vs3d = am%vs3d * (1 + gradient_s)
      am%vp3d = empirical_vp(am%vs3d)
      am%rho3d = empirical_rho(am%vp3d)
      call write_tmp_model()
    endif
    call synchronize_all()
    call sync_from_main_rank(am%vs3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%vp3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%rho3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
  end subroutine steepest_descent

  subroutine line_search(this)
    class (att_tomo), intent(inout) :: this
    real(kind=dp) :: chi, chi_local
    integer :: sit
    character(len=MAX_STRING_LEN) :: secname

    call write_log('Optimizing using halving stepping...',1,this%module)
    call write_gradient()
    do sit = 1, ap%inversion%max_sub_niter
      call write_log('Starting line search.',1,this%module)
      chi = 0
      call prepare_fwd_linesearch()
      do itype = 1, 2
        if (.not. ap%data%vel_type(itype)) cycle
        call select_type()
        call aq%forward_simulate(chi_local, .false., .false., verbose=.false.)
        chi = chi + chi_local
      enddo
      if (this%misfits(iter) > chi) then
        write(this%message, '(a,F0.4,a,F0.4)') 'Misfit reduced from ',this%misfits(iter),' to ',chi
        call write_log(this%message,0,this%module)
        write(this%message, '(a,F0.4," is ok, break sub-iterations for line search")') 'Step length of ',updatemax
        call write_log(this%message,1,this%module)
        exit
      else
        write(this%message, '(a,F0.4,a,F0.4)') 'Misfit increase from ',this%misfits(iter),' to ',chi
        call write_log(this%message,0,this%module)
        updatemax = updatemax * ap%inversion%maxshrink
        write(this%message, '(a,I2.2,a,F0.4)') 'In ',sit,'th sub-iteration, shrink step length to ',updatemax
        call write_log(this%message,1,this%module)
      endif
    enddo
    if (myrank == 0) then
      ! update & write model
      am%vs3d = am%vs3d_opt
      am%vp3d = empirical_vp(am%vs3d)
      am%rho3d = empirical_rho(am%vp3d)
      call write_tmp_model()
    endif
    call synchronize_all()
    call sync_from_main_rank(am%vs3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%vp3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%rho3d, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
  end subroutine line_search

  subroutine prepare_fwd_linesearch()
    ! class(att_tomo), intent(inout) :: this
    real(kind=dp) :: max_gk
    real(kind=dp), dimension(:,:,:), allocatable :: gradient_ls

    if (myrank == 0) then
      max_gk = maxval(abs(gradient_s))
      gradient_ls = -updatemax * gradient_s / max_gk
      am%vs3d_opt = am%vs3d * (1 + gradient_ls)
      am%vp3d_opt = empirical_vp(am%vs3d_opt)
      am%rho3d_opt = empirical_rho(am%vp3d_opt)    
    endif
    call synchronize_all()
    call sync_from_main_rank(am%vs3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%vp3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%rho3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
  end subroutine prepare_fwd_linesearch

  subroutine write_tmp_model()
    character(len=MAX_STRING_LEN) :: secname

    if (ap%output%verbose_level > 0) then
      write(secname,'(a,i3.3)') '/vs_',iter
      call h%add(secname, transpose_3(am%vs3d))
    endif

  end subroutine write_tmp_model

  subroutine write_gradient()
    character(len=MAX_STRING_LEN) :: secname

    if (ap%output%verbose_level > 0) then
      write(secname,'(a,i3.3)') '/gradient_',iter-1  
      call h%add(secname, transpose_3(gradient_s))
      do itype = 1, 2
        if (.not. ap%data%vel_type(itype)) cycle
        call select_type()
        write(secname,'(a,a,"_",i3.3)') '/kdensity_',trim(ap%data%gr_name(itype)),iter-1
        call h%add(secname, transpose_3(aq%ker_density))
      enddo
    endif

  end subroutine write_gradient

  subroutine select_type()
    if (itype == 1) then
      aq => aq_ph
      ap%data%igr = 0
    else
      aq => aq_gr
      ap%data%igr = 1
    endif
  end subroutine select_type

  subroutine break_iter(this, isbreak)
    class(att_tomo), intent(inout) :: this
    real(kind=dp) :: misfit_prev, misfit_curr, misfit_diff
    logical, intent(out) :: isbreak

    isbreak = .false.
    if (myrank == 0 .and.  iter > iter_store) then
      misfit_prev = sum(this%misfits(iter-iter_store:iter-1))/iter_store
      misfit_curr = sum(this%misfits(iter-iter_store+1:iter))/iter_store
      misfit_diff = (misfit_prev-misfit_curr)/misfit_prev
      write(this%message,'(a,i2,a,F6.4)') 'Misfit change of last', &
            iter_store,' iterations: ',misfit_diff
      call write_log(this%message,1,this%module)
      if (abs(misfit_diff) < ap%inversion%min_derr) isbreak = .true.
    endif
    call synchronize_all()
    call bcast_all(isbreak)
    
  end subroutine break_iter

  subroutine reset_opt(this)
    class(att_tomo), intent(inout) :: this
  
    if (myrank == 0) then
      if (.not. this%is_fwd) gradient_s = zeros(am%n_xyz(1),am%n_xyz(2),am%n_xyz(3))
      am%vs3d_opt = am%vs3d
      am%vp3d_opt = am%vp3d
      am%rho3d_opt = am%rho3d
    endif
    call synchronize_all()
    call sync_from_main_rank(am%vs3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%vp3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    call sync_from_main_rank(am%rho3d_opt, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
  end subroutine reset_opt

end module
