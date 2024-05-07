!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: Apl 2023, Initialize Codes
!
!=====================================================================

program rotate_model
  use hdf5_interface
  use csv_module
  use utils
  use constants
  use sph2loc
  use argparse, only: argparse_rotate_model

  implicit none
  
  character(len=MAX_STRING_LEN) :: input_file, output_file
  real(kind=dp) :: angle, center(2)
  real(kind=dp), dimension(:), allocatable :: lon, lat, dep, new_lat, new_lon
  real(kind=dp), dimension(:,:,:), allocatable :: vs3d
  real(kind=dp), dimension(:,:), allocatable :: vs_xyz
  type(hdf5_file) :: hf
  type(csv_file) :: fp
  logical :: status_ok
  integer :: i, j, k, nx, ny, nz

  call argparse_rotate_model(input_file, angle, center, output_file)

  ! Read the input file
  call hf%open(input_file)
  call hf%get('/lon', lon)
  call hf%get('/lat', lat)
  call hf%get('/dep', dep)
  call hf%get('/vs', vs3d)
  call hf%close()
  vs3d = transpose_3(vs3d)

  ! flatten vs3d
  nx = size(lon)
  ny = size(lat)
  nz = size(dep)
  allocate(vs_xyz(nx*ny*nz, 4))
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        vs_xyz((i-1)*ny*nz + (j-1)*nz + k, 1) = lat(j)
        vs_xyz((i-1)*ny*nz + (j-1)*nz + k, 2) = lon(i)
        vs_xyz((i-1)*ny*nz + (j-1)*nz + k, 3) = dep(k)
        vs_xyz((i-1)*ny*nz + (j-1)*nz + k, 4) = vs3d(i,j,k)
      end do
    end do
  end do

  ! rotate the model
  if (.not. (isnan(center(1)) .or. isnan(center(2)))) then
    call rtp_rotation_reverse(vs_xyz(:, 1), vs_xyz(:, 2), center(1), center(2), -angle, new_lat, new_lon)
    vs_xyz(:, 1) = new_lat
    vs_xyz(:, 2) = new_lon
  endif

  ! write to csv
  call fp%initialize(enclose_strings_in_quotes=.false.)
  call fp%open(output_file,n_cols=4,status_ok=status_ok)
  if ( .not. status_ok) then
    write(*,*) 'Cannot open '//trim(output_file)
    stop 
  endif
  call fp%add(['lat', 'lon', 'dep', 'vs '], trim_str=.true.)
  call fp%next_row()
  do i = 1, nx*ny*nz
    call fp%add(vs_xyz(i, 1), real_fmt='(F9.4)')
    call fp%add(vs_xyz(i, 2), real_fmt='(F9.4)')
    call fp%add(vs_xyz(i, 3), real_fmt='(F9.4)')
    call fp%add(vs_xyz(i, 4), real_fmt='(F9.4)')
    call fp%next_row()
  end do
  call fp%close(status_ok)

end program rotate_model
