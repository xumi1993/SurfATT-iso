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
! These codes come from Prof. Tongping's traveltime adjoint tomography program
! They are converted into subroutine, in order to be called by SurfATT
!
!=====================================================================
module sph2loc
  use constants
  implicit none

  interface rtp2xyz
    module procedure rtp2xyz_1d, rtp2xyz_sc
  end interface rtp2xyz

  interface rtp_rotation_reverse
    module procedure rtp_rotation_reverse_1, rtp_rotation_reverse_2
  end interface rtp_rotation_reverse

contains
  subroutine rtp2xyz_1d(r, theta, phi, x, y, z)
    real(kind=dp), dimension(:), intent(in) :: theta, phi
    real(kind=dp), intent(in) :: r
    real(kind=dp), dimension(:), allocatable, intent(out) :: x, y, z
    integer :: n

    n = size(theta)
    allocate(x(n),y(n),z(n))
    x = r * cosd(theta) * cosd(phi)
    y = r * cosd(theta) * sind(phi)
    z = r * sind(theta)
  end subroutine rtp2xyz_1d

  subroutine rtp2xyz_sc(r, theta, phi, x, y, z)
    real(kind=dp), intent(in) :: theta, phi, r
    real(kind=dp), intent(out) :: x, y, z

    x = r * cosd(theta) * cosd(phi)
    y = r * cosd(theta) * sind(phi)
    z = r * sind(theta)
  end subroutine rtp2xyz_sc

  ! Cartesian coordinates to Spherical coordinate
  subroutine xyz2rtp(x,y,z, r,theta,phi)
    ! theta: -90~90;  phi: -180~180
    real(kind=dp), dimension(:), intent(in) :: x, y, z
    real(kind=dp), dimension(:), allocatable, intent(out) :: r, theta, phi
    integer :: n, i

    n = size(x)
    allocate(r(n),theta(n),phi(n))
    r = sqrt(x**2+y**2+z**2)
    theta = asind(z/r)
    phi = asind(y/r/cosd(theta))
    do i = 1, n
      if(phi(i) > 0 .and. x(i)*y(i) < 0) phi(i) = 180 - phi(i)
      if(phi(i) < 0 .and. x(i)*y(i) > 0) phi(i) = -180 - phi(i)
    enddo
  end subroutine xyz2rtp    

  ! anti-clockwise rotation along x-axis
  subroutine rotate_x(x,y,z,theta)
    real(kind=dp), dimension(:), intent(inout) :: x, y, z
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(:), allocatable :: new_x, new_y, new_z
    integer :: n

    n = size(x)
    allocate(new_x(n),new_y(n),new_z(n))
    new_y = y *  cosd(theta) + z * (-sind(theta))
    new_z = y *  sind(theta) + z *  cosd(theta)
    y = new_y
    z = new_z
  end subroutine rotate_x
      
  ! anti-clockwise rotation along y-axis
  subroutine rotate_y(x,y,z,theta)
    real(kind=dp), dimension(:), intent(inout) :: x, y, z
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(:), allocatable :: new_x, new_y, new_z
    integer :: n

    n = size(x)
    allocate(new_x(n),new_y(n),new_z(n))
    new_x = x *  cosd(theta) + z * sind(theta)
    new_z = x * (-sind(theta)) + z * cosd(theta)
    x = new_x
    z = new_z
  end subroutine rotate_y

  ! anti-clockwise rotation along z-axis
  subroutine rotate_z(x,y,z,theta)
    real(kind=dp), dimension(:), intent(inout) :: x, y, z
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(:), allocatable :: new_x, new_y, new_z
    integer :: n

    n = size(x)
    allocate(new_x(n),new_y(n),new_z(n))
    new_x = x *  cosd(theta) + y * (-sind(theta))
    new_y = x *  sind(theta) + y *  cosd(theta)
    x = new_x
    y = new_y
  end subroutine rotate_z

  ! rotate to the new coordinate, satisfying the center r0,t0,p0 -> r0,0,0 and a anticlockwise angle psi
  subroutine rtp_rotation(t,p,theta0,phi0,psi,new_t, new_p)
    real(kind=dp), dimension(:), intent(in) :: t, p
    real(kind=dp), intent(in) :: theta0, phi0, psi
    real(kind=dp), dimension(:), allocatable, intent(out) :: new_t, new_p
    real(kind=dp), dimension(:), allocatable :: x, y, z, new_r
    integer :: n

    n = size(t)
    ! step 1: r,t,p -> x,y,z
    call rtp2xyz(1.0_dp,t,p,x,y,z)

    ! step 2: anti-clockwise rotation with -phi0 along z-axis:   r0,t0,p0 -> r0,t0,0
    call rotate_z(x,y,z,-phi0)

    ! step 3: anti-clockwise rotation with theta0 along y-axis:  r0,t0,0 -> r0,0,0
    call rotate_y(x,y,z,theta0)

    ! step 4: anti-clockwise rotation with psi along x-axis
    call rotate_x(x,y,z,psi)

    ! step 5: x,y,z -> r,t,p
    call xyz2rtp(x,y,z, new_r,new_t,new_p)
    
  end subroutine rtp_rotation

  subroutine rtp_rotation_reverse_1(new_t,new_p,theta0,phi0,psi,t,p)
    real(kind=dp), dimension(:), intent(in) :: new_t, new_p
    real(kind=dp), intent(in) :: theta0, phi0, psi
    real(kind=dp), dimension(:), allocatable, intent(out) :: t, p
    real(kind=dp), dimension(:), allocatable :: x, y, z, new_r
    integer :: n

    n = size(new_t)
    ! step 1: r,t,p -> x,y,z
    call rtp2xyz(1.0_dp,new_t,new_p,x,y,z)

    ! step 2: anti-clockwise rotation with -psi along x-axis
    call rotate_x(x,y,z,-psi)

    ! step 3: anti-clockwise rotation with -theta0 along y-axis:  r0,0,0 -> r0,t0,0 
    call rotate_y(x,y,z,-theta0)

    ! step 4: anti-clockwise rotation with phi0 along z-axis:   r0,t0,0 -> r0,t0,p0 
    call rotate_z(x,y,z,phi0)

    ! step 5: x,y,z -> r,t,p
    call xyz2rtp(x,y,z,new_r,t,p)
  end subroutine rtp_rotation_reverse_1
    
  subroutine rtp_rotation_reverse_2(new_t,new_p,theta0,phi0,psi,t,p)
    real(kind=dp), dimension(:,:), intent(in) :: new_t, new_p
    real(kind=dp), intent(in) :: theta0, phi0, psi
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: t, p
    real(kind=dp), dimension(:), allocatable :: x, y, z, new_r,new_t_flat,&
                                                new_p_flat, t_flat, p_flat
    integer :: nx, ny, i, j

    nx = size(new_t,1)
    ny = size(new_t,2)
    
    new_t_flat = reshape(new_t,[nx*ny])
    new_p_flat = reshape(new_p,[nx*ny])
    
    call rtp_rotation_reverse_1(new_t_flat,new_p_flat,theta0,phi0,psi,t_flat,p_flat)
    
    t = reshape(t_flat,[nx,ny])
    p = reshape(p_flat,[nx,ny])
  end subroutine rtp_rotation_reverse_2

end module sph2loc