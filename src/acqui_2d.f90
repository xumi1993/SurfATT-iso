module acqui_2d
  use utils
  use para
  use grid, ag_ph => att_grid_global_ph, ag_gr => att_grid_global_gr
  use src_rec, sr_ph => src_rec_global_ph, sr_gr => src_rec_global_gr
  ! use stdlib_math, only: linspace

  type, public :: att_acqui_2d
    real(kind=dp), dimension(:,:,:), pointer               :: svel, adj_s, &
                                                              ker_s, m11,m12,m22
    type(att_grid), pointer                                :: ag
    type(SrcRec), pointer                                  :: sr
    character(len=MAX_STRING_LEN)                          :: model_fname, module='ACQUI2D',&
                                                              final_fname, gr_name
    integer                                                :: nsrc, istart, iend, itype, iter = 0
    real(kind=dp)                                          :: updatemax, chi0
    real(kind=dp), dimension(:), pointer                   :: misfits
    integer, dimension(:,:), allocatable                   :: isrcs
    contains
    procedure :: init => att_acqui_2d_init, add_pert => att_acqui_2d_add_pert, &
                 write_model => att_acqui_2d_write_model, write_iter => att_acqui_2d_write_iter, &
                  write_obj_func => att_acqui_2d_write_obj_func
    procedure :: prepare_fwd, init_model, prepare_inv, scatter_src_gather, prepare_fwd_linesearch
    procedure, private :: construct_1d_ref_model, allocate_shm_arrays
  end type

  type(att_acqui_2d), target                               :: acqui_ph, acqui_gr
  integer :: win_adj_s, win_misfit, win_svel_acqui, win_xi, win_eta
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
    this%model_fname = trim(ap%output%output_path)//"/"//trim(modfile)
    this%final_fname = trim(ap%output%output_path)//"/"//trim('final_model.h5')
    
    call this%allocate_shm_arrays()
    call this%init_model()
  end subroutine att_acqui_2d_init

  subroutine scatter_src_gather(this)
    class(att_acqui_2d), intent(inout) :: this
    integer, dimension(:,:), allocatable :: isrcs

    isrcs = zeros(this%sr%npath, 2)
    this%nsrc = 0
    do i = 1, this%sr%nperiod
      do j = 1, this%sr%stations%nsta
        if (any(this%sr%periods_all==this%sr%periods(i)) .and. &
            any(this%sr%evtname==this%sr%stations%staname(j))) then          
          this%nsrc = this%nsrc+1
          isrcs(this%nsrc, 1) = i
          isrcs(this%nsrc, 2) = j
        endif
      enddo
    enddo
    this%isrcs = zeros(this%nsrc, 2)
    this%isrcs(1:this%nsrc, :) = isrcs(1:this%nsrc, :)
    call scatter_all_i(this%nsrc,mysize,myrank,this%istart,this%iend)
    call synchronize_all()
  end subroutine scatter_src_gather

  subroutine att_acqui_2d_add_pert(this,nx,ny,pert_vel,hmarg)
    class(att_acqui_2d), intent(inout) :: this
    integer, intent(in) :: nx,ny
    real(kind=dp), intent(in), optional :: pert_vel,hmarg
    real(kind=dp) :: ashmarg,aspert_vel,amp
    real(kind=dp), dimension(:), allocatable :: x_pert, y_pert
    integer :: i, j, ntaperx, ntapery

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
  end subroutine att_acqui_2d_add_pert

  subroutine construct_1d_ref_model(this)
    class(att_acqui_2d), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: zgrids, vs1d,svel
    integer :: i, nz

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
    if (myrank == 0) this%ag%svel = this%svel
    call synchronize_all()
  end subroutine prepare_fwd

  subroutine prepare_fwd_linesearch(this, xi_new, eta_new)
    class(att_acqui_2d), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), allocatable :: update_s, update_xi, update_eta
    real(kind=dp), dimension(:,:,:), allocatable :: xi_new, eta_new
    
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
    if (myrank == 0) then
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
    if (myrank == 0) this%adj_s = 0.0_dp
    call synchronize_all()
  end subroutine prepare_inv

  subroutine att_acqui_2d_write_iter(this)
    class(att_acqui_2d), intent(inout) :: this
    character(MAX_NAME_LEN) :: gr_name, secname

    if (myrank == 0) then
      gr_name = trim(ap%data%gr_name(itype))
      if(this%iter == 0) then
        call h5write(this%model_fname, '/stx_'//trim(gr_name), this%sr%stations%stx)
        call h5write(this%model_fname, '/sty_'//trim(gr_name), this%sr%stations%sty)
        call h5write(this%model_fname, '/periods_'//trim(gr_name), this%ag%periods)
        call h5write(this%model_fname, '/x', this%ag%xgrids)
        call h5write(this%model_fname, '/y', this%ag%ygrids)
      endif
      write(secname,'(a,i3.3)') '/vel_'//trim(gr_name)//'_',this%iter 
      call h5write(this%model_fname, secname, this%svel)
    endif
    call synchronize_all()
  end subroutine att_acqui_2d_write_iter

  subroutine att_acqui_2d_write_model(this)
    class(att_acqui_2d), intent(inout) :: this
    character(MAX_NAME_LEN) :: gr_name, secname

    if (myrank == 0) then
      call h5write(this%final_fname, '/stx_'//trim(this%gr_name), this%sr%stations%stx)
      call h5write(this%final_fname, '/sty_'//trim(this%gr_name), this%sr%stations%sty)
      call h5write(this%final_fname, '/periods_'//trim(this%gr_name), this%ag%periods)
      call h5write(this%final_fname, '/x', this%ag%xgrids)
      call h5write(this%final_fname, '/y', this%ag%ygrids)
      call h5write(this%final_fname, '/vel_'//trim(this%gr_name), this%svel)
    endif
    call synchronize_all()
  end subroutine att_acqui_2d_write_model

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
    call prepare_shm_array_dp_3d(this%ker_s, this%ag%nperiod, this%ag%nx, this%ag%ny, win_adj_s)
  end subroutine allocate_shm_arrays

end module acqui_2d