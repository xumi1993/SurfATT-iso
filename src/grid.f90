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
module grid
  use para, ap => att_para_global
  use src_rec, sr_gr => src_rec_global_gr, sr_ph => src_rec_global_ph
  use model, am => att_model_global
  use decomposer, amd => att_mesh_decomposer_global
  use topo
  use utils
  use surfker, only: fwdsurf3d_decomp
  use setup_att_log

  implicit none

  type, public :: att_grid
    real(kind=dp), dimension(:,:,:), pointer :: a, b, c,topo_angle
    real(kind=dp), dimension(:,:,:), pointer :: m11, m12, m22, ref_t, svel
    real(kind=dp), dimension(:), pointer :: periods, xgrids, ygrids
    real(kind=dp) :: dx, dy
    integer :: nperiod, igr, nx, ny
    character(len=MAX_STRING_LEN) :: module = "GRID"
    contains
    procedure :: init => init_grid, get_topo, fwdsurf
  end type

  integer :: win_topo, win_a, win_b, win_c, win_ref_t, win_svel, win_periods, &
             win_xgrids, win_ygrids, win_m11, win_m12, win_m22, win_topo_angle
  type(att_grid), target, public :: att_grid_global_ph, att_grid_global_gr

  contains

  subroutine init_grid(this, itype)
    class(att_grid), intent(inout) :: this
    integer, optional, intent(in) :: itype
    type(srcrec), pointer :: sr
    integer :: i
    
    if (present(itype)) then
      if (itype==1) then
        sr => sr_ph
        this%igr = 0
      elseif(itype==2) then
        sr => sr_gr
        this%igr = 1
      endif
    else
      sr => sr_ph
      this%igr = 0
    endif
    this%nperiod = sr%nperiod
    this%nx = am%n_xyz(1)
    this%ny = am%n_xyz(2)
    this%dx = am%d_xyz(1)
    this%dy = am%d_xyz(2)
    call prepare_shm_array_dp_1d(this%periods, this%nperiod, win_periods)
    call prepare_shm_array_dp_1d(this%xgrids, this%nx, win_xgrids)
    call prepare_shm_array_dp_1d(this%ygrids, this%ny, win_ygrids)
    call prepare_shm_array_dp_3d(this%a, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_a)
    call prepare_shm_array_dp_3d(this%b, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_b)
    call prepare_shm_array_dp_3d(this%c, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_c)
    call prepare_shm_array_dp_3d(this%topo_angle, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_topo_angle)
    call prepare_shm_array_dp_3d(this%m11, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_m11)
    call prepare_shm_array_dp_3d(this%m12, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_m12)
    call prepare_shm_array_dp_3d(this%m22, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_m22)
    call prepare_shm_array_dp_3d(this%ref_t, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_ref_t)
    call prepare_shm_array_dp_3d(this%svel, sr%nperiod, am%n_xyz(1), am%n_xyz(2), win_svel)

    if (local_rank == 0) then
      this%periods = sr%periods
      this%xgrids = am%xgrids
      this%ygrids = am%ygrids
      this%ref_t = 1._dp
    endif
    call synchronize_all()
  end subroutine init_grid

  subroutine get_topo(this)
    class(att_grid), intent(inout) :: this
    type(att_topo) :: at
    real(kind=dp), dimension(:), allocatable :: tmp, periods
    real(kind=dp), dimension(:,:), allocatable :: tmpto, fx, fy
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_a, tmp_b, tmp_c, tmp_agl
    real(kind=dp) :: sigma
    character(len=MAX_STRING_LEN) :: name
    integer :: win_topo, igr, ip, ix, iy, istart, iend

    if (ap%topo%is_consider_topo) then
      call write_log("Reading topography file",1, this%module)
      call at%read_parallel(ap%topo%topo_file)
      call at%grid_topo(am%xgrids, am%ygrids)
      allocate(tmp(this%nperiod))
      tmp_a = zeros(this%nperiod, this%nx, this%ny)
      tmp_b = zeros(this%nperiod, this%nx, this%ny)
      tmp_c = zeros(this%nperiod, this%nx, this%ny)
      tmp_agl = zeros(this%nperiod, this%nx, this%ny)
      call scatter_all_i(this%nperiod, mysize, myrank, istart, iend)
      if (myrank == 0) then
        call fwdsurf1d(am%vs1d,ap%data%iwave,&
                     igr,this%periods,&
                     am%zgrids,tmp)
      endif
      call synchronize_all()
      call bcast_all(tmp, this%nperiod)
      if (iend - istart >= 0) then
        do ip = istart, iend
          sigma = tmp(ip) * this%periods(ip) * ap%topo%wavelen_factor*km2deg
          tmpto = at%smooth(sigma)
          tmp_agl(ip,:,:) = at%calc_dip_angle(tmpto)
          call gradient_2_geo(tmpto, am%xgrids, am%ygrids, fx, fy)
          tmp_a(ip,:,:) = (1+fy**2) / (1 + fx**2 + fy**2)
          tmp_b(ip,:,:) = (1+fx**2) / (1 + fx**2 + fy**2)
          tmp_c(ip,:,:) = fx*fy / (1 + fx**2 + fy**2)
        enddo
      endif
      call synchronize_all()
      call sum_all(tmp_a, this%a, this%nperiod, this%nx, this%ny)
      call sum_all(tmp_b, this%b, this%nperiod, this%nx, this%ny)
      call sum_all(tmp_c, this%c, this%nperiod, this%nx, this%ny)
      call sum_all(tmp_agl, this%topo_angle, this%nperiod, this%nx, this%ny)
      call sync_from_main_rank(this%a, this%nperiod, this%nx, this%ny)
      call sync_from_main_rank(this%b, this%nperiod, this%nx, this%ny)
      call sync_from_main_rank(this%c, this%nperiod, this%nx, this%ny)
      call sync_from_main_rank(this%topo_angle, this%nperiod, this%nx, this%ny)
    else
      if (local_rank == 0) then
        this%a = 1.
        this%b = 1.
        this%c = 0.
      endif
    endif
    if (local_rank == 0) then
      this%m11 = this%a
      this%m22 = this%b
      this%m12 = -this%c
    endif
    call synchronize_all()

  end subroutine get_topo

  subroutine fwdsurf(this, vs3d)
    class(att_grid), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: vs3d
    real(kind=dp), dimension(:,:,:), allocatable :: tmpvel
    ! calculate surfave wave velocity

    call fwdsurf3d_decomp(vs3d,amd%loc_ix_start,amd%loc_ix_end,amd%loc_iy_start,amd%loc_iy_end,&
                          ap%data%iwave,this%igr,this%periods,am%zgrids,tmpvel) 
    call synchronize_all()
    call amd%collect_grid(tmpvel, this%svel)
    call sync_from_main_rank(this%svel, this%nperiod, this%nx, this%ny)
    call synchronize_all()
    
  end subroutine fwdsurf
end module