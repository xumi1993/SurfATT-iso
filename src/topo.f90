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
module topo
  use shared_par
  use hdf5_interface
  use utils
  use sph2loc
  use my_mpi

  implicit none

  type, public :: att_topo
    real(kind=dp), dimension(:), allocatable :: lon, lat
    real(kind=dp), dimension(:,:), allocatable :: z
    integer, dimension(2) :: dims
    real(kind=dp) :: dx, dy
    contains
    procedure :: read => read_topo, smooth => gaussian_smooth, &
                 read_parallel => read_topo_mpi, &
                 grid_topo, calc_dip_angle, rotate, write
  end type att_topo
  
  integer :: win_to_lat, win_to_lon, win_to_z
  
  contains

  subroutine read_topo_mpi(this, fname)
    class(att_topo), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: fname
    type(hdf5_file) :: h
    
    if (myrank == 0) then
      call this%read(fname)
    end if
    call synchronize_all()
    call bcast_all(this%dims, 2)
    if (myrank /= 0) then
      allocate(this%lon(this%dims(1)), this%lat(this%dims(2)),&
               this%z(this%dims(1), this%dims(2)))
    end if
    call bcast_all(this%lon, this%dims(1))
    call bcast_all(this%lat, this%dims(2))
    call bcast_all(this%z, this%dims(1)*this%dims(2))
    call bcast_all(this%dx)
    call bcast_all(this%dy)
  end subroutine read_topo_mpi

  subroutine read_topo(this, fname)
    class(att_topo), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: fname
    type(hdf5_file) :: h
    
    call h%open(fname, status='old', action='r')
    call h%get('/lon', this%lon)
    call h%get('/lat', this%lat)
    call h%get('/z', this%z)
    this%z = this%z/1000
    this%dx = this%lon(2) - this%lon(1)
    this%dy = this%lat(2) - this%lat(1)
    call h%close()
    this%dims = [size(this%lon), size(this%lat)]
  end subroutine read_topo

  function gaussian_smooth(this, sigma) result(topo)
    class(att_topo), intent(inout) :: this
    real(kind=dp), intent(in) :: sigma ! in deg
    real(kind=dp), dimension(:,:), allocatable :: topo

    topo = gaussian_smooth_geo_2(this%z, this%lon, this%lat, sigma)

  end function gaussian_smooth    

  subroutine grid_topo(this, xgrids, ygrids)
    class(att_topo), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: xgrids, ygrids
    real(kind=dp), dimension(:,:), allocatable :: xx, yy, topo_tmp
    integer :: i, j, k, ierr
    real(kind=dp) :: dist, az, baz, delta, sigma3, w, s
    integer :: nx, ny, count
    
    nx = size(xgrids)
    ny = size(ygrids)
    this%dx = xgrids(2) - xgrids(1)
    this%dy = ygrids(2) - ygrids(1)
    if (myrank == 0) then
      if(xgrids(1) < this%lon(1) .or. xgrids(size(xgrids)) > this%lon(this%dims(1)) .or. &
        ygrids(1) < this%lat(1) .or. ygrids(size(ygrids)) > this%lat(this%dims(2)) ) then
        print *, 'Error: Grid points are out of range'
        stop
      end if
      call meshgrid2(xgrids, ygrids, xx, yy)
      topo_tmp = transpose(interp2(this%lat, this%lon, this%z, yy, xx))
      deallocate(this%z, this%lat, this%lon)
      this%dims = [nx, ny]
      this%lon = xgrids
      this%lat = ygrids
      this%z = topo_tmp
    else
      deallocate(this%z, this%lat, this%lon)
      allocate(this%lon(nx), this%lat(ny), this%z(nx, ny))  
    endif
    call synchronize_all()
    call bcast_all(this%dims, 2)
    call bcast_all(this%lon, this%dims(1))
    call bcast_all(this%lat, this%dims(2))
    call bcast_all(this%z, this%dims(1)*this%dims(2))
  end subroutine grid_topo

  function calc_dip_angle(this, topo) result(angle)
    class(att_topo), intent(inout) :: this
    real(kind=dp), dimension(:,:), intent(in) :: topo
    real(kind=dp), dimension(:,:), allocatable :: angle
    real(kind=dp), dimension(:,:),allocatable :: Tx, Ty

    call gradient_2_geo(topo, this%lon, this%lat, Tx, Ty)
    angle = atand(sqrt(Tx**2 + Ty**2))
  end function calc_dip_angle
  
  subroutine rotate(this, xmin, xmax, ymin, ymax, clat, clon, angle)
    class(att_topo), intent(inout) :: this
    real(kind=dp), intent(in) :: angle, clat, clon, xmin, xmax, ymin, ymax
    real(kind=dp), dimension(:), allocatable :: x, y
    real(kind=dp), dimension(:,:), allocatable :: xx, yy, topo_tmp, xx_bk, yy_bk
    integer :: nx, ny

    x = arange(xmin, xmax, this%dx)
    y = arange(ymin, ymax, this%dy)
    call meshgrid2(x, y, xx, yy)
    call rtp_rotation_reverse(yy, xx, clat, clon, angle, yy_bk, xx_bk)
    if (minval(xx_bk) < this%lon(1) .or. maxval(xx_bk) > this%lon(this%dims(1)) .or. &
        minval(yy_bk) < this%lat(1) .or. maxval(yy_bk) > this%lat(this%dims(2))) then
      print *, 'Error: Grid points are out of range'
      print *, 'Rotated: ', minval(xx_bk), maxval(xx_bk), minval(yy_bk), maxval(yy_bk)
      print *, 'Original: ',this%lon(1), this%lon(this%dims(1)), this%lat(1), this%lat(this%dims(2))
      stop
    end if
    topo_tmp = interp2(this%lat, this%lon, this%z, yy_bk, xx_bk)
    deallocate(this%z, this%lat, this%lon)
    this%dims = [size(x), size(y)]
    this%lon = x
    this%lat = y
    this%z = transpose(topo_tmp)
  end subroutine rotate

  subroutine write(this, fname)
    class(att_topo), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: fname
    type(hdf5_file) :: h

    call h%open(fname, status='replace')
    call h%add('/lon', this%lon)
    call h%add('/lat', this%lat)
    call h%add('/z', this%z)
    call h%close()

  end subroutine write

end module
