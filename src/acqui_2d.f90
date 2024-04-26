!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: Mar 2024, Initialize Codes
!
!=====================================================================
module acqui_2d
  use utils
  use para
  use hdf5_interface
  use grid, ag_ph => att_grid_global_ph, ag_gr => att_grid_global_gr
  use src_rec, sr_ph => src_rec_global_ph, sr_gr => src_rec_global_gr
  ! use stdlib_math, only: linspace
  
  implicit none

  type, public :: att_acqui_2d
    real(kind=dp), dimension(:,:,:), pointer               :: svel, adj_s,&
                                                              ker_s, m11,m12,m22
    real(kind=dp), dimension(:,:,:), allocatable           :: adj_s_local
    type(att_grid), pointer                                :: ag
    type(SrcRec), pointer                                  :: sr
    character(len=MAX_STRING_LEN)                          :: model_fname, module='ACQUI2D',&
                                                              final_fname, gr_name
    integer                                                :: itype, iter = 0
    real(kind=dp)                                          :: updatemax, chi0
    real(kind=dp), dimension(:), pointer                   :: misfits
    type(hdf5_file)                                        :: h
    contains
    procedure :: init => att_acqui_2d_init, add_pert => att_acqui_2d_add_pert, &
                 write_model => att_acqui_2d_write_model, write_iter => att_acqui_2d_write_iter, &
                  write_obj_func => att_acqui_2d_write_obj_func, write_target_model => att_acqui_2d_write_target_model
    procedure :: prepare_fwd, init_model, prepare_inv, prepare_fwd_linesearch
    procedure, private :: construct_1d_ref_model, allocate_shm_arrays
  end type

  type(att_acqui_2d), target                               :: acqui_ph, acqui_gr
  integer :: win_adj_s, win_misfit, win_svel_acqui, win_isrcs_2d, win_ker_s
  character(len=MAX_STRING_LEN), private                   :: message

  contains

  subroutine att_acqui_2d_init(this, itype)
    class(att_acqui_2d), intent(inout) :: this
    integer, intent(in) :: itype

    this%itype = itype
    this%gr_name = trim(ap%data%gr_name(itype))
    if (itype == 1) then
      this%ag => ag_ph
      this%sr => sr_ph
    elseif(itype == 2) then
      this%ag => ag_gr
      this%sr => sr_gr
    endif
    this%model_fname = trim(ap%output%output_path)//"/"//'model_iter_'//trim(this%gr_name)//'.h5'
    this%final_fname = trim(ap%output%output_path)//"/"//'final_model_'//trim(this%gr_name)//'.h5'
    
    call this%allocate_shm_arrays()
    call this%init_model()
  end subroutine att_acqui_2d_init

  ! subroutine scatter_src_gather(this)
  !   class(att_acqui_2d), intent(inout) :: this
  !   integer, dimension(:,:), allocatable :: isrcs
  !   integer, dimension(:), allocatable :: iperiods
  !   integer :: i, j, np

  !   isrcs = zeros(this%sr%npath, 2)
  !   iperiods = zeros(this%sr%nperiod)
  !   this%nsrc = 0
  !   if (local_rank == 0) then
  !     do j = 1, this%sr%stations%nsta
  !       if (any(this%sr%evtname==this%sr%stations%staname(j))) then  
  !         call this%sr%get_periods_by_src(this%sr%stations%staname(j), iperiods, np)
  !         do i = 1, np
  !           this%nsrc = this%nsrc+1        
  !           isrcs(this%nsrc, 1) = iperiods(i)
  !           isrcs(this%nsrc, 2) = j
  !         enddo
  !       endif
  !     enddo
  !     write(message,'(a,i0," ",a,a,i0,a)') 'Scatter ',this%nsrc,&
  !           trim(this%gr_name),' events to ',mysize," processors"
  !     call write_log(message,1,this%module)
  !   endif
  !   call synchronize_all()
  !   call bcast_all(this%nsrc)
  !   call prepare_shm_array_i_2d(this%isrcs, this%nsrc, 2, win_isrcs_2d)
  !   if (local_rank == 0) this%isrcs(1:this%nsrc, :) = isrcs(1:this%nsrc, :)
  !   call synchronize_all()
  !   call scatter_all_i(this%nsrc,mysize,myrank,this%istart,this%iend)
  ! end subroutine scatter_src_gather

  subroutine att_acqui_2d_add_pert(this,nx,ny,pert_vel,hmarg)
    class(att_acqui_2d), intent(inout) :: this
    integer, intent(in) :: nx,ny
    real(kind=dp), intent(in), optional :: pert_vel,hmarg
    real(kind=dp) :: ashmarg,aspert_vel,amp
    real(kind=dp), dimension(:), allocatable :: x_pert, y_pert
    integer :: i, j, k, ntaperx, ntapery

    if (present(pert_vel)) then 
      aspert_vel = pert_vel 
    else
      aspert_vel = 0.08
    endif

    if (present(hmarg)) then
      ashmarg = hmarg
    else
      ashmarg = 0.
    endif
    
    call write_log('Adding perturbation to the model', 1, this%module)
    write(message, '(a,2i3)') 'Perturbation grid num: ', nx, ny
    call write_log(message, 1, this%module)
    write(message, '(a,f6.3,a,f6.3,a,f6.3)') 'Perturbation of velocity: ', aspert_vel
    call write_log(message, 1, this%module)
    if (myrank == 0) then
      call this%construct_1d_ref_model()
      do i = 1, this%ag%nperiod
        ntaperx = ashmarg/this%ag%dx
        ntapery = ashmarg/this%ag%dy
        x_pert = zeros(this%ag%nx)
        x_pert(ntaperx+1:this%ag%ny-ntaperx) = &
          sin(nx*pi*arange(this%ag%nx-2*ntaperx)/(this%ag%nx-2*ntaperx))
        y_pert = zeros(this%ag%ny)
        y_pert(ntapery+1:this%ag%ny-ntapery) = &
          sin(ny*pi*arange(this%ag%ny-2*ntapery)/(this%ag%ny-2*ntapery))
        do j = 1, this%ag%nx
          do k = 1, this%ag%ny
            amp = aspert_vel * x_pert(j) * y_pert(k)
            this%svel(i,j,k) = this%svel(i,j,k) * (1 + amp)
          enddo
        enddo
      enddo
    endif
    call synchronize_all()
    call sync_from_main_rank(this%svel, this%ag%nperiod, this%ag%nx, this%ag%ny)
  end subroutine att_acqui_2d_add_pert

  subroutine construct_1d_ref_model(this)
    class(att_acqui_2d), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: zgrids, vs1d,svel
    integer :: i, j, nz

    zgrids = arange(ap%domain%depth(1),ap%domain%depth(2),ap%domain%interval(3))
    nz = size(zgrids)
    vs1d = linspace(ap%inversion%vel_range(1), ap%inversion%vel_range(2), nz)
    svel = zeros(this%sr%nperiod)
    call fwdsurf1d(vs1d,ap%data%iwave,ap%data%igr,&
                   this%sr%periods,zgrids,svel)
    do i = 1, this%ag%nx
      do j = 1, this%ag%ny
        this%svel(:,i,j) = svel
      enddo
    enddo
  end subroutine construct_1d_ref_model

  subroutine prepare_fwd(this)
    class(att_acqui_2d), intent(inout) :: this
    if (local_rank == 0) this%ag%svel = this%svel
    call synchronize_all()
  end subroutine prepare_fwd

  subroutine prepare_fwd_linesearch(this)
    class(att_acqui_2d), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), allocatable :: update_s
    
    if (myrank == 0) then
      update_s = zeros(this%ag%nperiod, this%ag%nx, this%ag%ny)
      update_s = this%updatemax*this%ker_s/maxval(abs(this%ker_s))
      this%ag%svel = this%svel*(1+update_s)
    endif
    call synchronize_all()
  end subroutine prepare_fwd_linesearch

  subroutine init_model(this)
    class(att_acqui_2d), intent(inout) :: this
    integer :: i, j
    if (local_rank == 0) then
      do i = 1, this%ag%nx
        do j = 1, this%ag%ny
          this%svel(:,i,j) = this%sr%meanvel
        enddo
      enddo
      this%ag%svel = this%svel
    endif
    call synchronize_all()
  end subroutine init_model

  subroutine prepare_inv(this)
    class(att_acqui_2d), intent(inout) :: this
    this%adj_s_local = zeros(this%ag%nperiod, this%ag%nx, this%ag%ny)
    if (local_rank == 0) then
      this%ker_s = 0._dp
      this%adj_s = 0._dp
      this%sr%tt_fwd = 0._dp
    endif
    call synchronize_all()
  end subroutine prepare_inv

  subroutine att_acqui_2d_write_iter(this)
    class(att_acqui_2d), intent(inout) :: this
    character(MAX_NAME_LEN) :: gr_name, secname

    if (myrank == 0) then
      gr_name = trim(ap%data%gr_name(this%itype))
      if(this%iter == 0) then
        call this%h%add('/stlo', this%sr%stations%stlo)
        call this%h%add('/stla', this%sr%stations%stla)
        call this%h%add('/period', this%ag%periods)
        call this%h%add('/lon', this%ag%xgrids)
        call this%h%add('/lat', this%ag%ygrids)
      endif
      write(secname,'(a,i3.3)') '/vel'//trim(gr_name)//'_',this%iter 
      call this%h%add(secname, transpose_3(this%svel))
    endif
  end subroutine att_acqui_2d_write_iter

  subroutine att_acqui_2d_write_model(this)
    class(att_acqui_2d), intent(inout) :: this
    character(MAX_NAME_LEN) :: gr_name, secname
    type(hdf5_file) :: hf

    if (myrank == 0) then
      call hf%open(this%final_fname, status='new', action='write')

      call hf%add('/stlo', this%sr%stations%stlo)
      call hf%add('/stla', this%sr%stations%stla)
      call hf%add('/period', this%ag%periods)
      call hf%add('/lon', this%ag%xgrids)
      call hf%add('/lat', this%ag%ygrids)
      call hf%add('/vel', transpose_3(this%svel))
      call hf%close()
    endif
    call synchronize_all()
  end subroutine att_acqui_2d_write_model

  subroutine att_acqui_2d_write_target_model(this)
    class(att_acqui_2d), intent(inout) :: this
    character(MAX_NAME_LEN) :: gr_name, secname
    character(len=MAX_STRING_LEN) :: target_model_fname
    type(hdf5_file) :: hf

    if (myrank == 0) then
      target_model_fname =  trim(ap%output%output_path)//"/target_model_"//trim(this%gr_name)//".h5"
      call hf%open(target_model_fname, status='new', action='write')

      call hf%add('/period', this%ag%periods)
      call hf%add('/lon', this%ag%xgrids)
      call hf%add('/lat', this%ag%ygrids)
      call hf%add('/vel', transpose_3(this%svel))
      call hf%close()
    endif
    call synchronize_all()
  end subroutine att_acqui_2d_write_target_model

  subroutine att_acqui_2d_write_obj_func(this)
    class(att_acqui_2d), intent(inout) :: this

    if (myrank == 0) then
      ! write objective function
      write(IOUT, '(i3,",",f10.2,",",f6.4)') &
            this%iter-1,this%misfits(this%iter), this%updatemax
      call flush(IOUT)
    endif
    call synchronize_all()
  end subroutine att_acqui_2d_write_obj_func

  subroutine allocate_shm_arrays(this)
    class(att_acqui_2d), intent(inout) :: this
    call prepare_shm_array_dp_1d(this%misfits, ap%inversion%niter, win_misfit)
    call prepare_shm_array_dp_3d(this%svel, this%sr%nperiod, this%ag%nx, this%ag%ny, win_svel_acqui)
    call prepare_shm_array_dp_3d(this%ker_s, this%ag%nperiod, this%ag%nx, this%ag%ny, win_ker_s)
    call prepare_shm_array_dp_3d(this%adj_s, this%ag%nperiod, this%ag%nx, this%ag%ny, win_adj_s)
  end subroutine allocate_shm_arrays

end module acqui_2d