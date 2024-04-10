!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                 (there are currently many more authors!)
!                           (c) October 2023
!   
!     Changing History: Oct 2023, Initialize Codes
!
!=====================================================================

module model
  use my_mpi
  use para, ap => att_para_global
  use src_rec, staall => stations_global, sr_ph=>src_rec_global_ph, sr_gr=>src_rec_global_gr
  use decomposer, amd => att_mesh_decomposer_global
  use utils
  use hdf5_interface
  use surfker, only: fwdsurf1d, depthkernel1d
  use setup_att_log
  implicit none

  integer :: win_vs3d, win_vp3d, win_rho3d,win_vs1d, &
             win_vs3d_opt, win_vp3d_opt, win_rho3d_opt
  real(kind=dp), dimension(:), private, allocatable              :: anch, n_pi
  
  type, public :: att_model
    real(kind=dp)                                          :: d_xyz(3)
    integer                                                :: n_xyz(3)
    real(kind=dp), dimension(:), allocatable               :: xgrids, ygrids, zgrids
    real(kind=dp), dimension(:), pointer                   :: vs1d
    real(kind=dp), dimension(:,:), allocatable             :: xinv, yinv, zinv, topo
    real(kind=dp), dimension(:,:,:), pointer               :: vs3d, vp3d, rho3d
    real(kind=dp), dimension(:,:,:), pointer               :: vs3d_opt, vp3d_opt, rho3d_opt
    integer, dimension(:,:), allocatable                   :: igrid
    character(MAX_NAME_LEN), private                       :: module='MODEL'
    ! integer                                                :: grid_istart, grid_iend
    ! integer, dimension(:), allocatable, public             :: grid_local_index
    contains
    procedure :: init => initialize_model, write => write_model
    procedure :: inv1d, get_init_model, get_inv_grids, add_pert, prepare_model_opt
  end type
  
  type(att_model), public                                  :: att_model_global

  contains

  subroutine initialize_model(this)
    class(att_model), intent(inout) :: this
    real(kind=dp) :: xbeg, xend, ybeg, yend
    integer :: i, j, n
    character(len=MAX_STRING_LEN) :: msg    

    this%d_xyz = ap%domain%interval
    this%zgrids = arange(ap%domain%depth(1),ap%domain%depth(2),this%d_xyz(3))
    xbeg = minval(staall%stlo) - ap%domain%num_grid_margin*this%d_xyz(1)
    xend = maxval(staall%stlo) + ap%domain%num_grid_margin*this%d_xyz(1)
    ybeg = minval(staall%stla) - ap%domain%num_grid_margin*this%d_xyz(2)
    yend = maxval(staall%stla) + ap%domain%num_grid_margin*this%d_xyz(2)
    this%xgrids = arange(xbeg, xend, this%d_xyz(1))
    this%ygrids = arange(ybeg, yend, this%d_xyz(2))
    this%n_xyz = [size(this%xgrids), size(this%ygrids), size(this%zgrids)]
    write(msg, '(a,3i4)') 'Model grids: nx,ny,nz: ', this%n_xyz
    call write_log(msg,1,this%module)
    write(msg, '(a,f0.4,", ",f0.4)') 'Lon range: ', xbeg, xend
    call write_log(msg,1,this%module)
    write(msg, '(a,f0.4,", ",f0.4)') 'Lat range: ', ybeg, yend
    call write_log(msg,1,this%module)
    call this%get_inv_grids()
    call amd%init(this%n_xyz(1), this%n_xyz(2))
    call synchronize_all()

  end subroutine initialize_model

  subroutine get_init_model(this)
    class(att_model), intent(inout) :: this
    integer :: i, j, ier, niter
    character(len=:), allocatable :: msger
    real(kind=dp), dimension(:), allocatable :: misfits
    real(kind=dp), dimension(:,:,:), allocatable :: vstmp

    call prepare_shm_array_dp_3d(this%vs3d, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_vs3d)
    call prepare_shm_array_dp_3d(this%vp3d, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_vp3d)
    call prepare_shm_array_dp_3d(this%rho3d, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_rho3d)
    call prepare_shm_array_dp_1d(this%vs1d, this%n_xyz(3), win_vs1d)
    ! this%vs1d = linspace(ap%inversion%vel_range(1), ap%inversion%vel_range(2), this%n_xyz(3))
    if (myrank == 0) then
      if (ap%inversion%init_model_type == 0) then
        this%vs1d = linspace(ap%inversion%vel_range(1), ap%inversion%vel_range(2), this%n_xyz(3))
      elseif (ap%inversion%init_model_type == 1) then
        call this%inv1d(ap%inversion%niter, this%vs1d, niter, misfits)
      elseif (ap%inversion%init_model_type == 2) then
        ! call load_npy(ap%inversion%init_model_path, vstmp, ier, msger)
        call h5read(ap%inversion%init_model_path, '/vs', vstmp)
        if (any(shape(vstmp) /= this%n_xyz)) then
          write(*,*) 'Shape of '//trim(ap%inversion%init_model_path)//' dose not match with' //&
                     ' shape of computational domain.'
          stop
        endif
        do i = 1, this%n_xyz(3)
          this%vs1d(i) = sum(vstmp(:,:,i))/(this%n_xyz(1)*this%n_xyz(2))
        enddo
        this%vs3d = vstmp
      endif
      if (ap%inversion%init_model_type < 2) then
        do i = 1, this%n_xyz(1)
          do j = 1, this%n_xyz(2)
            this%vs3d(i, j, :) = this%vs1d
          enddo
        enddo
      endif
      this%vp3d = empirical_vp(this%vs3d)
      this%rho3d = empirical_rho(this%vp3d)
    endif
    call synchronize_all()
  end subroutine get_init_model

  subroutine get_inv_grids(this)
    class(att_model), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: ixgrids, iygrids, izgrids, zref,zadd
    integer :: i, ninvx, ninvy, ninvz, nset, polar

    ninvx = ap%inversion%n_inv_grid(1)
    ninvy = ap%inversion%n_inv_grid(2)
    ninvz = ap%inversion%n_inv_grid(3)
    nset = ap%inversion%ncomponents
    
    ixgrids = linspace( &
      this%xgrids(1) - (this%xgrids(this%n_xyz(1))-this%xgrids(1))/(ninvx-2),&
      this%xgrids(this%n_xyz(1)) + (this%xgrids(this%n_xyz(1))-this%xgrids(1))/(ninvx-2),&
      ninvx)
    ! ixgrids = linspace( &
    !   this%xgrids(1) - (this%xgrids(this%n_xyz(1))-this%xgrids(1))/(ninvx-2),&
    !   this%xgrids(this%n_xyz(1)), ninvx)
    iygrids = linspace( &
      this%ygrids(1) - (this%ygrids(this%n_xyz(2))-this%ygrids(1))/(ninvy-2),&
      this%ygrids(this%n_xyz(2)), ninvy)
    izgrids = linspace( &
      this%zgrids(1) - 0.5*(this%zgrids(this%n_xyz(3))-this%zgrids(1))/(ninvz-2), &
      this%zgrids(this%n_xyz(3)) + 0.5*(this%zgrids(this%n_xyz(3))-this%zgrids(1))/(ninvz-2),&
      ninvz)
    zref = zeros(ninvz)
    zadd = zeros(ninvz)
    zref(1:ninvz-1) = izgrids(2:)
    zref(ninvz) = izgrids(ninvz)
    zadd = (zref-izgrids)/(nset+1)
    zadd(1) = 0
    this%xinv = zeros(ninvx, nset)
    this%yinv = zeros(ninvy, nset)
    this%zinv = zeros(ninvz, nset)
    polar = 1
    do i = 1, nset
      polar = -polar
      this%xinv(:, i) = ixgrids+polar*(i-1)*(ixgrids(2)-ixgrids(1))/(nset+1)
      ! this%xinv(:, i) = ixgrids+i*(ixgrids(2)-ixgrids(1))/(nset+1)
      this%yinv(:, i) = iygrids+(i-1)*(iygrids(2)-iygrids(1))/(nset+1)
      this%zinv(:, i) = izgrids+(i-1)*zadd
    enddo
  end subroutine get_inv_grids

  subroutine inv1d(this, iter_num, vsinv, niter, misfits)
    class(att_model), intent(inout) :: this
    type(srcrec), pointer :: sr
    integer, intent(in) :: iter_num
    real(kind=dp), parameter :: minderr = 0.0001
    real(kind=dp), dimension(this%n_xyz(3)), intent(out) :: vsinv
    integer, intent(out) :: niter
    real(kind=dp), dimension(:), allocatable, intent(out) :: misfits
    real(kind=dp), dimension(this%n_xyz(3)) :: update, sen,update_total
    real(kind=dp), dimension(:,:), allocatable :: sen_vs, sen_vp, sen_rho
    character(len=MAX_STRING_LEN) :: msg
    real(kind=dp) :: derr, chi, sigma
    integer :: iter, ip, itype
    real(kind=dp), dimension(:),allocatable :: tmp

    vsinv = this%vs1d
    misfits = zeros(iter_num)
    sigma = this%zgrids(this%n_xyz(3))/ap%inversion%n_inv_grid(3)/2
    call write_log('Do 1D inverison using averaged surface wave data',1,this%module)
    do iter = 1, iter_num
      update_total = 0.
      do itype = 1, 2
        if (.not. ap%data%vel_type(itype)) cycle
        if (itype == 1) then
          sr => sr_ph
          ap%data%igr = 0
        else
          sr => sr_gr
          ap%data%igr = 1
        endif
        sen_vs = zeros(sr%nperiod,this%n_xyz(3))
        sen_vp = zeros(sr%nperiod,this%n_xyz(3))
        sen_rho = zeros(sr%nperiod,this%n_xyz(3))
        tmp = zeros(sr%nperiod)
        call fwdsurf1d(vsinv,ap%data%iwave,ap%data%igr,&
                       sr%periods,this%zgrids,tmp)
        chi = 0.5*sum((sr%meanvel-tmp)**2)
        misfits(iter) = misfits(iter) + chi
        call depthkernel1d(vsinv,this%n_xyz(3),ap%data%iwave,&
                            ap%data%igr,sr%nperiod,&
                            sr%periods,this%zgrids,&
                            sen_vs, sen_vp, sen_rho)
        update = 0.
        do ip = 1, sr%nperiod
          sen = sen_vs(ip,:) + sen_vp(ip,:)*dalpha_dbeta(vsinv) + &
                sen_rho(ip,:)*drho_dalpha(empirical_vp(vsinv))*dalpha_dbeta(vsinv)
          update = update - sen * (sr%meanvel(ip)-tmp(ip))
        enddo
        update = update / sr%nperiod
        update = smooth_1(update, this%zgrids, sigma)
        update_total = update_total + update*ap%data%weights(itype)
      enddo
      write(msg,'(a,i0,a,f8.4)') 'Iteration ', iter, ' misfit: ', misfits(iter)
      call write_log(msg,0,this%module)
      if (iter > 1) then
        derr = abs(misfits(iter)-misfits(iter-1))
        if (derr < minderr) exit
      endif
      vsinv = vsinv * (1-update_total)
    enddo
    niter = iter
    
  end subroutine inv1d

  subroutine add_pert(this,nx,ny,nz,pert_vel,hmarg,anom_size)
    class(att_model), intent(inout) :: this
    integer,  intent(in) :: nx, ny, nz
    real(kind=dp), intent(in), optional :: pert_vel,hmarg,anom_size
    real(kind=dp) :: ashmarg,aspert_vel,amp, asanom_size
    real(kind=dp), dimension(:), allocatable :: x_pert, y_pert, z_pert
    real(kind=dp), dimension(:,:,:), allocatable :: vs_pert
    real(kind=dp), dimension(3) :: para
    integer :: ntaperx, ntapery, i, j, k

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

    if (present(anom_size)) then
      asanom_size = anom_size
    else
      asanom_size = 0.
    endif

    if (myrank == 0) then
      ntaperx = ashmarg/this%d_xyz(1)
      ntapery = ashmarg/this%d_xyz(2)
      x_pert = zeros(this%n_xyz(1))
      x_pert(ntaperx+1:this%n_xyz(1)-ntaperx) = &
        sin(nx*pi*arange(this%n_xyz(1)-2*ntaperx)/(this%n_xyz(1)-2*ntaperx))
      y_pert = zeros(this%n_xyz(2))
      y_pert(ntapery+1:this%n_xyz(2)-ntapery) = &
        sin(ny*pi*arange(this%n_xyz(2)-2*ntapery)/(this%n_xyz(2)-2*ntapery))
      z_pert = zeros(this%n_xyz(3))
      if (asanom_size /= 0) then
        para = dep_anom(this,nz,asanom_size)
        z_pert = sin(2*pi*((sqrt(para(1)**2+para(2)*arange(0,this%n_xyz(3)-1))-para(1))/para(3)))
      else
        z_pert = sin(nz*pi*arange(this%n_xyz(3))/this%n_xyz(3))
      endif      
      vs_pert = zeros(this%n_xyz(1), this%n_xyz(2), this%n_xyz(3))
      do i = 1, this%n_xyz(1)
        do j = 1, this%n_xyz(2)
          do k = 1, this%n_xyz(3)
            vs_pert(i,j,k) = x_pert(i)*y_pert(j)*z_pert(k)*aspert_vel
          enddo
        enddo
      enddo
      this%vs3d = this%vs3d*(1+vs_pert)
      this%vp3d = empirical_vp(this%vs3d)
      this%rho3d = empirical_rho(this%vp3d)
    endif
    call synchronize_all()  
  end subroutine add_pert

  subroutine write_model(this, subname)
    class(att_model), intent(inout) :: this
    real(kind=dp),dimension(:,:), allocatable :: xyzdat 
    character(len=*), intent(in) :: subname
    character(len=MAX_STRING_LEN) :: fname
    type(hdf5_file) :: h
    integer :: i
    
    if (myrank == 0) then
      fname = trim(ap%output%output_path)//&
                    '/'//trim(subname)//'.h5'
      call h%open(fname, status='new', action='write')
      call h%add('/x',this%xgrids)
      call h%add('/y',this%ygrids)
      call h%add('/z',this%zgrids)
      call h%add('/vs',this%vs3d)
      call h%close()
    endif
    call synchronize_all()    
  end subroutine write_model

  subroutine dep_anom_fun(maxanchor,npara,para,fitfun)
    integer :: maxanchor,npara,i
    real(kind=dp), dimension(3) :: para
    real(kind=dp) :: fitfun(maxanchor),anomfun(maxanchor)
    
    do i = 1,maxanchor
        anomfun(i) = (sqrt(para(1)**2+para(2)*(anch(i)-1))-para(1))/para(3)
    end do
    fitfun(:) = abs(anomfun(:) - n_pi(1:maxanchor))
  end subroutine dep_anom_fun
  
  function dep_anom(this,nz,asanom_size) result(para)
    use minpack
    class(att_model), intent(inout) :: this
    integer :: nz, info=0, i
    real(kind=dp), dimension(3) :: para
    real(kind=dp) :: asanom_size, tol=1D-7, anom_size_inc, nanomtop
    real(kind=dp), dimension(:), allocatable :: fitfun
    character(len=MAX_STRING_LEN) :: msg

    nanomtop = asanom_size/(this%zgrids(2)-this%zgrids(1))
    anom_size_inc = (this%n_xyz(3)-1-nz*nanomtop)/(2*nz*(nz-1))
    fitfun = zeros(nz*2+1)
    anch = zeros(nz*2+1)
    n_pi = zeros(nz*2+1)
    anch(1) = 1
    n_pi(1) = 0
    para = [1,1,1]
    write(msg,'(a,f6.2,a,f3.1,a)') 'Depth anomaly anchor: ', (anch(1)-1)*(this%zgrids(2)-this%zgrids(1)), 'km, ', n_pi(1)*2, 'pi'
    call write_log(trim(msg),0,this%module)
    do i=2,nz*2+1
      anch(i) = anch(i-1) + (nanomtop-anom_size_inc)/2 + max(i-2,0)*anom_size_inc
      n_pi(i) = n_pi(i-1) + 0.25
      write(msg,'(a,f6.2,a,f3.1,a)') 'Depth anomaly anchor: ', (anch(i)-1)*(this%zgrids(2)-this%zgrids(1)), 'km, ', n_pi(i)*2, 'pi'
      call write_log(trim(msg),0,this%module)
    end do
    call lmdif1(dep_anom_fun,nz*2+1,3,para(:),fitfun(1:nz*2+1),tol,info)
  end function dep_anom

  subroutine prepare_model_opt(this)
    class(att_model), intent(inout) :: this
    
    call prepare_shm_array_dp_3d(this%vs3d_opt, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_vs3d_opt)
    call prepare_shm_array_dp_3d(this%vp3d_opt, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_vp3d_opt)
    call prepare_shm_array_dp_3d(this%rho3d_opt, this%n_xyz(1), this%n_xyz(2), this%n_xyz(3), win_rho3d_opt)
    call synchronize_all()
  end subroutine prepare_model_opt

end module
