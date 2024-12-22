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
module utils
  use constants
  implicit none

  integer, public, parameter :: IPRE = 4
  integer, public, parameter :: RPRE = dp
  integer, public, parameter :: DPRE = dp

  interface zeros
    module procedure zeros1, zeros2, zeros3, zeros4
  end interface zeros

  interface ones
    module procedure ones1, ones2, ones3
  end interface ones
  public :: ones
  private :: ones1, ones2, ones3

  interface randu
    module procedure randu0, randu1, randu2, randu3
  end interface randu

  interface randn
    module procedure randn0, randn1, randn2, randn3
  end interface randn

  interface diff
    module procedure diff1, diff2
  end interface diff
  public :: diff
  private :: diff1, diff2
  
  interface linspace
    module procedure linspace_r8r8, linspace_i4i4, &
                     linspace_r8i4,  linspace_i4r8
  end interface linspace

  interface arange
    module procedure arange0, arange1_i4, arange1_r8, arange_i4, arange_r8
  end interface arange

  interface interp1
    module procedure interp1_0, interp1_1
  end interface interp1

  interface interp2
    module procedure interp2_0_dp, interp2_1_dp, interp2_2_dp
  end interface interp2

  interface interp3
    module procedure interp3_0_dp
  end interface

  interface append
    module procedure append_i, append_r8, append_ch_name
  end interface append

  interface find_loc
    module procedure find_str_loc, find_int_loc
  end interface find_loc

  interface meshgrid_ij
    module procedure meshgrid2_ij, meshgrid3_ij
  end interface meshgrid_ij

  interface gps2dist
    module procedure gps2dist_scalar, gps2dist_1, gps2dist_2
  end interface gps2dist

  interface empirical_vp
    module procedure empirical_vp_1, empirical_vp_2, empirical_vp_3
  end interface empirical_vp
  
  interface empirical_rho
    module procedure empirical_rho_1, empirical_rho_2, empirical_rho_3
  end interface

  interface drho_dalpha
    module procedure drho_dalpha_0, drho_dalpha_1, drho_dalpha_2, drho_dalpha_3
  end interface drho_dalpha

  interface dalpha_dbeta
    module procedure dalpha_dbeta_1, dalpha_dbeta_2, dalpha_dbeta_3
  end interface dalpha_dbeta

contains
function drho_dalpha_0(vp) result(rho)
    double precision, intent(in) :: vp
    double precision :: rho
    rho = 1.6612 - 2*0.4721*vp + 3*0.0671*vp**2 - 4*0.0043*vp**3 + 5*0.000106*vp**4
end function

function drho_dalpha_1(vp) result(rho)
    double precision, dimension(:), intent(in) :: vp
    double precision, dimension(:), allocatable :: rho
    allocate(rho(size(vp)))
    rho = 1.6612 - 2*0.4721*vp + 3*0.0671*vp**2 - 4*0.0043*vp**3 + 5*0.000106*vp**4
end function

function drho_dalpha_2(vp) result(rho)
    double precision, dimension(:,:), intent(in) :: vp
    double precision, dimension(:,:), allocatable :: rho
    allocate(rho(size(vp,1), size(vp,2)))
    rho = 1.6612 - 2*0.4721*vp + 3*0.0671*vp**2 - 4*0.0043*vp**3 + 5*0.000106*vp**4
end function

function drho_dalpha_3(vp) result(rho)
    double precision, dimension(:,:,:),intent(in) :: vp
    double precision, dimension(:,:,:), allocatable :: rho
    allocate(rho(size(vp,1), size(vp,2), size(vp,3)))
    rho = 1.6612 - 2*0.4721*vp + 3*0.0671*vp**2 - 4*0.0043*vp**3 + 5*0.000106*vp**4
end function

function dalpha_dbeta_1(vs) result(dvp)
    double precision, dimension(:),intent(in) :: vs
    double precision, dimension(:), allocatable :: dvp
    allocate(dvp(size(vs)))
    dvp = 2.0947 - 2*0.8206*vs + 3*0.2683*vs**2 - 4*0.0251*vs**3
end function

function dalpha_dbeta_2(vs) result(dvp)
    double precision, dimension(:,:),intent(in) :: vs
    double precision, dimension(:,:), allocatable :: dvp
    allocate(dvp(size(vs,1), size(vs,2)))
    dvp = 2.0947 - 2*0.8206*vs + 3*0.2683*vs**2 - 4*0.0251*vs**3
end function

function dalpha_dbeta_3(vs) result(dvp)
    double precision, dimension(:,:,:),intent(in) :: vs
    double precision, dimension(:,:,:), allocatable :: dvp
    allocate(dvp(size(vs,1), size(vs,2), size(vs,3)))
    dvp = 2.0947 - 2*0.8206*vs + 3*0.2683*vs**2 - 4*0.0251*vs**3
end function

function empirical_vp_1(vs) result(vp)
    double precision, dimension(:),intent(in) :: vs
    double precision, dimension(:), allocatable :: vp
    allocate(vp(size(vs)))
    vp = 0.9409 + 2.0947*vs - 0.8206*vs**2+ 0.2683*vs**3 - 0.0251*vs**4
end function

function empirical_vp_2(vs) result(vp)
    double precision, dimension(:,:),intent(in) :: vs
    double precision, dimension(:,:), allocatable :: vp
    allocate(vp(size(vs,1),size(vs,2)))
    vp = 0.9409 + 2.0947*vs - 0.8206*vs**2+ 0.2683*vs**3 - 0.0251*vs**4
end function

function empirical_vp_3(vs) result(vp)
    double precision, dimension(:,:,:),intent(in) :: vs
    double precision, dimension(:,:,:), allocatable :: vp
    allocate(vp(size(vs,1),size(vs,2),size(vs,3)))
    vp = 0.9409 + 2.0947*vs - 0.8206*vs**2+ 0.2683*vs**3 - 0.0251*vs**4
end function

function empirical_rho_1(vp) result(rho)
  double precision, dimension(:), intent(in) :: vp
  double precision, dimension(:), allocatable :: rho
  allocate(rho(size(vp)))
  rho=1.6612*vp - 0.4721*vp**2 + &
        0.0671*vp**3 - 0.0043*vp**4 + & 
        0.000106*vp**5
end function

function empirical_rho_2(vp) result(rho)
  double precision, dimension(:,:), intent(in) :: vp
  double precision, dimension(:,:), allocatable :: rho
  allocate(rho(size(vp,1), size(vp,2)))
  rho=1.6612*vp - 0.4721*vp**2 + &
        0.0671*vp**3 - 0.0043*vp**4 + & 
        0.000106*vp**5
end function

function empirical_rho_3(vp) result(rho)
  double precision, dimension(:,:,:), intent(in) :: vp
  double precision, dimension(:,:,:), allocatable :: rho
  allocate(rho(size(vp,1), size(vp,2), size(vp,3)))
  rho=1.6612*vp - 0.4721*vp**2 + &
        0.0671*vp**3 - 0.0043*vp**4 + & 
        0.000106*vp**5
end function

!=======================================================================
! zeros
!-----------------------------------------------------------------------
! zeros creates array all of zeros.
!
! Syntax
!-----------------------------------------------------------------------
! x = zeros(dim1)
! A = zeros(dim1, dim2)
! X = zeros(dim1, dim2, dim3)
!
! Description
!-----------------------------------------------------------------------
! x = zeros(dim1) returns a dim1 vector of zeros.
!
! A = zeros(dim1, dim2) returns a dim1-by-dim2 matrix of zeros.
!
! X = zeros(dim1, dim2, dim3) returns a dim1-by-dim2-by-dim3
! 3-dimensional matrix of zeros.
!
! Examples
!-----------------------------------------------------------------------
! x = zeros(3)
! x =
!     0.  0.  0.
!
! A = zeros(3, 3)
! A =
!     0.  0.  0.
!     0.  0.  0.
!     0.  0.  0.
!=======================================================================

  pure function zeros1(dim1)
    real(kind = RPRE), dimension(:), allocatable :: zeros1
    integer(kind = IPRE), intent(in) :: dim1
    integer(kind = IPRE) :: ierr

    allocate(zeros1(dim1), stat = ierr)
    ! if ( ierr .ne. 0 ) then
    !   stop "Error: in zeros, could not allocate array."
    ! else
    zeros1 = 0.0d0
    ! end if
    return
  end function zeros1

  pure function zeros2(dim1, dim2)
    real(kind = RPRE), dimension(:,:), allocatable :: zeros2
    integer(kind = IPRE), intent(in) :: dim1, dim2
    integer(kind = IPRE) :: ierr

    allocate(zeros2(dim1, dim2), stat = ierr)
    ! if ( ierr .ne. 0 ) then
    !   print *, "Error: in zeros, could not allocate array."
    !   stop
    ! else
    zeros2 = 0.0d0
    ! end if
    return
  end function zeros2

  pure function zeros3(dim1, dim2, dim3)
    real(kind = RPRE), dimension(:,:,:), allocatable :: zeros3
    integer(kind = IPRE), intent(in) :: dim1, dim2, dim3
    integer(kind = IPRE) :: ierr

    allocate(zeros3(dim1, dim2, dim3), stat = ierr)
    ! if ( ierr .ne. 0 ) then
    !   print *, "Error: in zeros, could not allocate array."
    !   stop
    ! else
    zeros3 = 0.0d0
    ! end if
    return
  end function zeros3

  pure function zeros4(dim1, dim2, dim3, dim4)
    real(kind = RPRE), dimension(:,:,:,:), allocatable :: zeros4
    integer(kind = IPRE), intent(in) :: dim1, dim2, dim3, dim4
    integer(kind = IPRE) :: ierr

    allocate(zeros4(dim1, dim2, dim3, dim4), stat = ierr)
    ! if ( ierr .ne. 0 ) then
    !   print *, "Error: in zeros, could not allocate array."
    !   stop
    ! else
    zeros4 = 0.0d0
    ! end if
    return
  end function zeros4

!=======================================================================
! ones
!-----------------------------------------------------------------------
! ones creates array all of ones.
!
! Syntax
!-----------------------------------------------------------------------
! x = ones(dim1)
! A = ones(dim1, dim2)
! X = ones(dim1, dim2, dim3)
!
! Description
!-----------------------------------------------------------------------
! x = ones(dim1) returns a dim1 vector of ones.
!
! A = ones(dim1, dim2) returns a dim1-by-dim2 matrix of ones.
!
! X = ones(dim1, dim2, dim3) returns a dim1-by-dim2-by-dim3
! 3-dimensional matrix of ones.
!
! Examples
!-----------------------------------------------------------------------
! x = ones(3)
! x =
!     1.  1.  1.
!
! A = ones(3, 3)
! A =
!     1.  1.  1.
!     1.  1.  1.
!     1.  1.  1.
!=======================================================================

  function ones1(dim1)
    real(kind = RPRE), dimension(:), allocatable :: ones1
    integer(kind = IPRE), intent(in) :: dim1
    
    allocate(ones1(dim1))
    ones1 = 1.0d0
    return
  end function ones1
  
  function ones2(dim1, dim2)
    real(kind = RPRE), dimension(:,:), allocatable :: ones2
    integer(kind = IPRE), intent(in) :: dim1, dim2
    
    allocate(ones2(dim1, dim2))
    ones2 = 1.0d0
    return
  end function ones2
  
  function ones3(dim1, dim2, dim3)
    real(kind = RPRE), dimension(:,:,:), allocatable :: ones3
    integer(kind = IPRE), intent(in) :: dim1, dim2, dim3
    
    allocate(ones3(dim1, dim2, dim3))
    ones3 = 1.0d0
    return
  end function ones3

!=======================================================================
! randu
!-----------------------------------------------------------------------
! randu generates uniformly distributed random numbers.
!
! Syntax
!-----------------------------------------------------------------------
! x = randu()
! x = randu(dim1)
! A = randu(dim1, dim2)
! X = randu(dim1, dim2, dim3)
!
! Description
!-----------------------------------------------------------------------
! x = randu() returns a single uniformly distributed random number in
! the interval [0,1].
!
! x = randu(dim1) returns a dim1 vector of uniformly distributed random
! numbers.
!
! A = randu(dim1, dim2) returns a dim1-by-dim2 matrix of uniformly
! distributed random numbers.
!
! X = randu(dim1, dim2, dim3) returns a dim1-by-dim2-by-dim3
! 3-dimensional matrix of uniformly distributed random numbers.
!
! Examples
!-----------------------------------------------------------------------
! x = randu()
!     0.383413825
!
! x = randu(5)*2 - 1
!     0.640258908  -0.873707294   0.787327528
!=======================================================================

  real(kind = DPRE) function randu0()
    call random_number(randu0)
    return
  end function randu0

  function randu1(dim1)
    real(kind = DPRE), dimension(:), allocatable :: randu1
    integer(kind = IPRE), intent(in) :: dim1

    allocate(randu1(dim1))
    call random_number(randu1)
    return
  end function randu1

  function randu2(dim1, dim2)
    real(kind = DPRE), dimension(:,:), allocatable :: randu2
    integer(kind = IPRE), intent(in) :: dim1, dim2

    allocate(randu2(dim1, dim2))
    call random_number(randu2)
    return
  end function randu2

  function randu3(dim1, dim2, dim3)
    real(kind = DPRE), dimension(:,:,:), allocatable :: randu3
    integer(kind = IPRE), intent(in) :: dim1, dim2, dim3

    allocate(randu3(dim1, dim2, dim3))
    call random_number(randu3)
    return
  end function randu3

!=======================================================================
! randn
!-----------------------------------------------------------------------
! randn generates normally distributed random numbers using polar
! Box-Muller algorithm.
!
! Syntax
!-----------------------------------------------------------------------
! x = randn()
! x = randn(dim1)
!
! Description
!-----------------------------------------------------------------------
! x = randn() returns a single normally distributed random number with
! mean 0 and standard deviation 1.
!
! x = randn(dim1) returns a dim1 vector of normally distributed random
! numbers.
!
! A = randn(dim1, dim2) returns a dim1-by-dim2 matrix of normally
! distributed random numbers.
!
! X = randn(dim1, dim2, dim3) returns a dim1-by-dim2-by-dim3
! 3-dimensional matrix of normally distributed random numbers.
!
! Examples
!-----------------------------------------------------------------------
! x = randn(3)
!     -1.22003853  -0.211721316   0.522971511
!=======================================================================

  real(kind = DPRE) function randn0()
    real(kind = RPRE) :: u, v, s

    do
      u = 2.*randu() - 1.
      v = 2.*randu() - 1.
      s = u*u + v*v
      if ( (s .gt. 0.) .and. (s .lt. 1.) ) exit
    end do
    randn0 = u * sqrt( -2.0d0 * log(s) / s )
    return
  end function randn0

  function randn1(dim1)
    real(kind = DPRE), dimension(:), allocatable :: randn1
    integer(kind = IPRE), intent(in) :: dim1
    integer(kind = IPRE) :: i

    allocate(randn1(dim1))
    do i = 1, dim1
      randn1(i) = randn0()
    end do
    return
  end function randn1

  function randn2(dim1, dim2)
    real(kind = DPRE), dimension(:,:), allocatable :: randn2
    integer(kind = IPRE), intent(in) :: dim1, dim2
    integer(kind = IPRE) :: i, j

    allocate(randn2(dim1, dim2))
    do i = 1, dim1
      do j = 1, dim2
        randn2(i,j) = randn()
      end do
    end do
    return
  end function randn2

  function randn3(dim1, dim2, dim3)
    real(kind = DPRE), dimension(:,:,:), allocatable :: randn3
    integer(kind = IPRE), intent(in) :: dim1, dim2, dim3
    integer(kind = IPRE) :: i, j, k

    allocate(randn3(dim1, dim2, dim3))
    do i = 1, dim1
      do j = 1, dim2
        do k = 1, dim3
          randn3(i,j,k) = randn()
        end do
      end do
    end do
    return
  end function randn3

!=======================================================================
! diff
!-----------------------------------------------------------------------
! diff computes differences of arrays
!
! Syntax
!-----------------------------------------------------------------------
! y = diff(x)
! y = diff(x, n)
! B = diff(A)
! B = diff(A, n)
! B = diff(A, dim)
! B = diff(A, n, dim)
!
! Description
!-----------------------------------------------------------------------
! y = diff(x) returns differences between adjacent elements of vector x.
!
! y = diff(x, n) returns the nth difference by applying the diff(x)
! operator recursively n times.
!
! B = diff(A) returns differences between adjacent elements of array A
! along the first dimension.
!
! B = diff(A, n) returns the nth difference by applying the diff(A)
! operator recursively n times.
!
! B = diff(A, dim) returns differences between adjacent elements of
! array A along the dimension given by dim.
!
! B = diff(A, n, dim) returns the nth difference along the dimension
! given by dim by applying the diff(A, dim) operator recursively
! n times.
!=======================================================================

  function diff1(x, n)
    real(kind = DPRE), dimension(:), allocatable :: diff1
    real(kind = DPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE), intent(in), optional :: n
    integer(kind = IPRE) :: opt_n, i

    opt_n = 1
    if (present(n)) opt_n = n

    diff1 = x
    do i = 1, opt_n
      diff1 = diff1(2:) - diff1(:size(diff1)-1)
    end do
    return
  end function diff1

  function diff2(A, n, dim)
    real(kind = DPRE), dimension(:,:), allocatable :: diff2
    real(kind = DPRE), dimension(:,:), intent(in) :: A
    integer(kind = IPRE), intent(in), optional :: n, dim
    integer(kind = IPRE) :: opt_n, i

    opt_n = 1
    if (present(n)) opt_n = n

    diff2 = A
    if ((.not. present(dim)) .or. (dim .eq. 1)) then
      do i = 1, opt_n
        diff2 = diff2(2:,:) - diff2(:size(diff2,1)-1,:)
      end do
    elseif (dim .eq. 2) then
      do i = 1, opt_n
        diff2 = diff2(:,2:) - diff2(:,:size(diff2,2)-1)
      end do
    end if
    return
  end function diff2
!=======================================================================
! linspace
!-----------------------------------------------------------------------
! linspace creates a linearly spaced vector.
!
! Syntax
!-----------------------------------------------------------------------
! x = linspace(x1, x2, n)
!
! Description
!-----------------------------------------------------------------------
! x = linspace(x1, x2, n) returns a vector of n evenly spaced points
! between x1 and x2.
!
! Examples
!-----------------------------------------------------------------------
! x = linspace(0, 10, 11)
!     0.  1.  2.  3.  4.  5.  6.  7.  8.  9.  10.
!=======================================================================

  pure function linspace_r8r8(first, last, n)
    real(kind = DPRE), dimension(:), allocatable :: linspace_r8r8
    real(kind = DPRE), intent(in) :: first, last
    integer(kind = IPRE), intent(in) :: n
    integer(kind = IPRE) :: i
    real(kind = DPRE) :: step

    allocate(linspace_r8r8(n))
    step = ( last - first ) / ( n-1 )
    linspace_r8r8 = first + step * real([ ( i-1, i = 1, n ) ], DPRE)
    linspace_r8r8(n) = last
    return
  end function linspace_r8r8

  pure function linspace_i4i4(first, last, n)
    real(kind = DPRE), dimension(:), allocatable :: linspace_i4i4
    integer, intent(in) :: first, last
    integer(kind = IPRE), intent(in) :: n

    linspace_i4i4 = linspace(real(first, kind = 8), real(last, kind = 8), n)
    return
  end function linspace_i4i4

  pure function linspace_r8i4(first, last, n)
    real(kind = DPRE), dimension(:), allocatable :: linspace_r8i4
    real(kind = DPRE), intent(in) :: first
    integer, intent(in) :: last
    integer(kind = IPRE), intent(in) :: n

    linspace_r8i4 = linspace(first, real(last, kind = DPRE), n)
    return
  end function linspace_r8i4

  pure function linspace_i4r8(first, last, n)
    real(kind = DPRE), dimension(:), allocatable :: linspace_i4r8
    integer, intent(in) :: first
    real(kind = DPRE), intent(in) :: last
    integer(kind = IPRE), intent(in) :: n

    linspace_i4r8 = linspace(real(first, kind = 8), last, n)
    return
  end function linspace_i4r8

  pure function arange_i4(first, last, step) result(x)
    real(kind = DPRE), dimension(:), allocatable :: x
    integer(kind = IPRE), intent(in) :: first, last, step
    integer(kind = IPRE) :: size, i

    size = int((last - first) / step) + 1
    x = [ ( first + i * step, i = 0, size - 1 ) ]
    return
  end function arange_i4

  pure function arange_r8(first, last, step) result(x)
    real(kind = DPRE), dimension(:), allocatable :: x
    real(kind = DPRE), intent(in) :: first, last, step
    integer(kind = IPRE) :: i, size
    
    size = int((last - first) / step) + 1
    x = [ ( first + i * step, i = 0, size - 1 ) ]
    return
  end function arange_r8

  pure function arange0(num) result(x)
    real(kind = DPRE), dimension(:), allocatable :: x
    integer(kind = IPRE), intent(in) :: num
    integer(kind = IPRE) :: i

    x = [ ( i, i = 1, num) ]
    return
  end function

  pure function arange1_i4(first, last) result(x)
    real(kind = DPRE), dimension(:), allocatable :: x
    integer(kind = IPRE), intent(in) :: first, last
    integer(kind = IPRE) :: size, i

    size = int((last - first)) + 1
    x = [ ( first + i, i = 0, size - 1 ) ]
    return
  end function

  pure function arange1_r8(first, last) result(x)
    real(kind = DPRE), dimension(:), allocatable :: x
    real(kind = DPRE), intent(in) :: first, last
    integer(kind = IPRE) :: size, i

    size = int((last - first)) + 1
    x = [ ( first + i, i = 0, size - 1 ) ]
    return
  end function

  function sort(x, order)
    real(kind = DPRE), dimension(:), allocatable :: sort
    real(kind = DPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE), intent(in), optional :: order
    integer(kind = IPRE) :: n

    n = size(x)
    sort = x
    if ((.not. present(order)) .or. (order .eq. 1)) then
      call quicksort(sort, n, 1)
    elseif (order .eq. 2) then
      call quicksort(sort, n, 2)
    end if
    return

  contains

    !-------------------------------------------------------------------
    ! quicksort
    !-------------------------------------------------------------------
    recursive subroutine quicksort(x, n, order)
      real(kind = DPRE), dimension(n), intent(inout) :: x
      integer(kind = IPRE), intent(in) :: n, order
      integer(kind = IPRE) :: left, right, marker
      real(kind = DPRE) :: pivot, tmp

      if (n .gt. 1) then
        left = 0
        right = n + 1
        pivot = x(randi0_0(n))

        select case(order)
          case(1)
            do while ( left .lt. right )
              left = left + 1
              right = right - 1
              do while ( x(left) .lt. pivot )
                left = left + 1
              end do
              do while ( x(right) .gt. pivot )
                right = right - 1
              end do
              if ( left .lt. right ) then
                tmp = x(left)
                x(left) = x(right)
                x(right) = tmp
              end if
            end do
          case(2)
            do while ( left .lt. right )
              left = left + 1
              right = right - 1
              do while ( x(left) .gt. pivot )
                left = left + 1
              end do
              do while ( x(right) .lt. pivot )
                right = right - 1
              end do
              if ( left .lt. right ) then
                tmp = x(left)
                x(left) = x(right)
                x(right) = tmp
              end if
            end do
        end select

        if ( left .eq. right ) then
          marker = left + 1
        else
          marker = left
        end if

        call quicksort(x(:marker-1), marker-1, order)
        call quicksort(x(marker:), n-marker+1, order)
      end if
      return
    end subroutine quicksort

  end function sort

  integer(kind = IPRE) function randi0_0(imax)
    integer(kind = IPRE), intent(in) :: imax

    randi0_0 = floor( randu0() * real(imax) ) + 1
    return
  end function randi0_0

!=======================================================================
! interp1
!-----------------------------------------------------------------------
! interp1 performs a linear interpolation.
!
! Syntax
!-----------------------------------------------------------------------
! vq = interp1(x, v, xq)
!
! Description
!-----------------------------------------------------------------------
! vq = interp1(x, v, xq) returns the evaluated vector yq at the query
! points in xq using a linear interpolation.
!=======================================================================

  function interp1_0(x, v, xq) result(vq)
    real(kind = DPRE) :: vq
    real(kind = DPRE), intent(in) :: xq
    real(kind = DPRE), dimension(:), intent(in) :: x, v
    integer(kind = IPRE) :: i, x1, x2, ix(2)
    real(kind = DPRE) :: vn, xr(2), vr(2)

    x1 = minloc(xq - x, 1, mask = xq .ge. x)
    x2 = maxloc(xq - x, 1, mask = xq .lt. x)
    if ( x2 .ne. 0 ) then
      vn = abs( (x(x2) - x(x1)) )
      xr = x( [ x1, x2 ] )
      vr = v( [ x1, x2 ] )
      vq = vr(1) * ( xr(2) - xq ) + vr(2) * ( xq - xr(1) )
      vq = vq / vn
    else
      vq = v(size(v))
    end if
    return
  end function interp1_0

  function interp1_1(x, v, xq) result(vq)
    real(kind = DPRE), dimension(:), allocatable :: vq
    real(kind = DPRE), dimension(:), intent(in) :: xq, x, v
    integer(kind = IPRE) :: i, n

    n = size(xq)
    vq = zeros(n)
    do i = 1, n
      vq(i) = interp1_0(x, v, xq(i))
    end do
    return
  end function interp1_1

  pure function transpose_3(x) result(y)
    real(kind = DPRE), dimension(:,:,:), intent(in) :: x
    real(kind = DPRE), dimension(size(x,3), size(x,2), size(x,1)) :: y
    integer(kind = IPRE) :: i, j, k

    do i = 1, size(x,1)
      do j = 1, size(x,2)
        do k = 1, size(x,3)
          y(k,j,i) = x(i,j,k)
        end do
      end do
    end do
    return
  end function transpose_3

  pure function transpose_4(x) result(y)
    real(kind = DPRE), dimension(:,:,:,:), intent(in) :: x
    real(kind = DPRE), dimension(size(x,4), size(x,3), size(x,2), size(x,1)) :: y
    integer(kind = IPRE) :: i, j, k, l

    do i = 1, size(x,1)
      do j = 1, size(x,2)
        do k = 1, size(x,3)
          do l = 1, size(x,4)
            y(l,k,j,i) = x(i,j,k,l)
          enddo
        end do
      end do
    end do
    return
  end function transpose_4

  pure subroutine meshgrid2_ij(ax, ay, x, y)
    real(kind = DPRE), dimension(:), intent(in) :: ax, ay
    real(kind = DPRE), dimension(:,:), allocatable, intent(out) :: x, y
    integer(kind = IPRE) :: i, m, n

    m = size(ax)
    n = size(ay)
    if (.not. allocated(x)) allocate(x(m, n))
    if (.not. allocated(y)) allocate(y(m, n))
    do i = 1, n
      x(:,i) = ax
    end do
    do i = 1, m
      y(i,:) = ay
    end do
    return
  end subroutine meshgrid2_ij

  pure subroutine meshgrid3_ij(ax, ay, az, x, y, z)
    real(kind = DPRE), dimension(:), intent(in) :: ax, ay, az
    real(kind = DPRE), dimension(:,:,:), allocatable, intent(out) :: x, y, z
    integer(kind = IPRE) :: i, m, n, j, k

    m = size(ax)
    n = size(ay)
    k = size(az)
    if (.not. allocated(x)) allocate(x(m, n, k))
    if (.not. allocated(y)) allocate(y(m, n, k))
    if (.not. allocated(z)) allocate(z(m, n, k))
    do i = 1, n
      do j = 1, k
        x(:,i,j) = ax
      enddo
    end do
    do i = 1, m
      do j = 1, k
        y(i,:,j) = ay
      enddo
    end do
    do i = 1, m
      do j = 1, n
        z(i, j, :) = az
      enddo
    enddo
    return
  end subroutine meshgrid3_ij

  pure subroutine meshgrid2(ax, ay, x, y)
    real(kind = DPRE), dimension(:), intent(in) :: ax, ay
    real(kind = DPRE), dimension(:,:), allocatable, intent(out) :: x, y
    integer(kind = IPRE) :: i, m, n

    m = size(ax)
    n = size(ay)
    if (.not. allocated(x)) allocate(x(n, m))
    if (.not. allocated(y)) allocate(y(n, m))
    do i = 1, n
      x(i,:) = ax
    end do
    do i = 1, m
      y(:,i) = ay
    end do
    return
  end subroutine meshgrid2

  pure function interp2_0_dp(x, y, v, xq, yq) result(vq)
    real(kind = DPRE) :: vq
    real(kind = DPRE), intent(in) :: xq, yq
    real(kind = DPRE), dimension(:), intent(in) :: x, y
    real(kind = DPRE), dimension(:,:), intent(in) :: v
    integer(kind = IPRE) :: i, x1, y1, x2, y2, ix(4), iy(4)
    real(kind = DPRE) :: vn, xr(2), yr(2), N(4), vr(4)

    x1 = minloc(xq - x, 1, mask = xq .ge. x)
    y1 = minloc(yq - y, 1, mask = yq .ge. y)
    x2 = maxloc(xq - x, 1, mask = xq .lt. x)
    y2 = maxloc(yq - y, 1, mask = yq .lt. y)
    vn = abs( (x(x2) - x(x1)) &
            * (y(y2) - y(y1)) )
    xr = x( [ x1, x2 ] )
    yr = y( [ y1, y2 ] )
    ix = [ 2, 1, 2, 1 ]
    iy = [ 2, 2, 1, 1 ]
    do i = 1, 4
      N(i) = abs( (xr(ix(i)) - xq) * (yr(iy(i)) - yq) )
    end do
    vr = reshape(v( [ x1, x2 ], &
                    [ y1, y2 ] ), shape = [ 4 ])
    vq = dot_product(vr, N/vn)
    return
  end function interp2_0_dp

  pure function interp2_1_dp(x, y, v, xq, yq) result(vq)
    real(kind = DPRE), dimension(:), allocatable :: vq
    real(kind = DPRE), dimension(:), intent(in) :: xq, yq, x, y
    real(kind = DPRE), dimension(:,:), intent(in) :: v
    integer(kind = IPRE) :: i, n

    n = size(xq)
    vq = zeros(n)
    do i = 1, n
      vq(i) = interp2_0_dp(x, y, v, xq(i), yq(i))
    end do
    return
  end function interp2_1_dp

  pure function interp2_2_dp(x, y, v, xq, yq) result(vq)
    real(kind = DPRE), dimension(:,:), allocatable :: vq
    real(kind = DPRE), dimension(:), intent(in) :: x, y
    real(kind = DPRE), dimension(:,:), intent(in) :: v, xq, yq
    integer(kind = IPRE) :: m, n

    m = size(xq, 1)
    n = size(xq, 2)
    vq = reshape( interp2_1_dp(y, x, v, [ yq ], [ xq ]), shape = [ m, n ] )
    return
  end function interp2_2_dp

  function interp3_0_dp(valx,valy,valz,valin,x,y,z) result(valout)


    !*** Lookup 1D tables with Cartesian coordinates
    real(kind=DPRE), dimension(:), intent(in) :: valx
    real(kind=DPRE), dimension(:), intent(in) :: valy
    real(kind=DPRE), dimension(:), intent(in) :: valz

    !*** Input grid of data
    real(kind=DPRE), dimension(:,:,:), intent(in) :: valin

    !*** Coordinate for interpolated value
    real(kind=DPRE), intent(in) :: x, y, z

    !*** Temporary variables
    real(kind=DPRE) :: xix1, xix2, yiy1, yiy2, ziz1, ziz2
    real(kind=DPRE) :: x2x1, x1x2, y2y1, y1y2, z2z1, z1z2
    real(kind=DPRE) :: facx1y1z1, facx1y1z2, facx1y2z1, facx1y2z2, facx2y1z1, facx2y1z2, facx2y2z1, facx2y2z2

    !*** Output value to be interpolated
    real(kind=DPRE) :: valout
    integer :: indx, indy, indz  !!! Previous guest
    integer :: kx, ky, kz, m, nx,ny,nz
    
    nx = size(valx)
    ny = size(valy)
    nz = size(valz)

    call locate_bissection(valx,nx,x,indx)
    call locate_bissection(valy,ny,y,indy)
    call locate_bissection(valz,nz,z,indz)

    m=2
    kx = min(max(indx-(m-1)/2,1),nx+1-m)
    ky = min(max(indy-(m-1)/2,1),ny+1-m)
    kz = min(max(indz-(m-1)/2,1),nz+1-m)


    !*** x_i - x1
    xix1 = x - valx(kx)
    yiy1 = y - valy(ky)
    ziz1 = z - valz(kz)

    !*** x_i - x2
    xix2 = x - valx(kx+1)
    yiy2 = y - valy(ky+1)
    ziz2 = z - valz(kz+1)

    !*** x1 - x2
    x1x2 = 1./ (valx(kx) - valx(kx+1))
    y1y2 = 1./ (valy(ky) - valy(ky+1))
    z1z2 = 1./ (valz(kz) - valz(kz+1))

    !*** x2 - x1
    x2x1 = 1./(valx(kx+1) - valx(kx))
    y2y1 = 1./(valy(ky+1) - valy(ky))
    z2z1 = 1./(valz(kz+1) - valz(kz))

    !*** Factors
    facx1y1z1 = xix2*yiy2*ziz2 * x1x2*y1y2*z1z2 * valin(kx,ky,kz)
    facx1y1z2 = xix2*yiy2*ziz1 * x1x2*y1y2*z2z1 * valin(kx,ky,kz+1)
    facx1y2z1 = xix2*yiy1*ziz2 * x1x2*y2y1*z1z2 * valin(kx,ky+1,kz)
    facx1y2z2 = xix2*yiy1*ziz1 * x1x2*y2y1*z2z1 * valin(kx,ky+1,kz+1)
    facx2y1z1 = xix1*yiy2*ziz2 * x2x1*y1y2*z1z2 * valin(kx+1,ky,kz)
    facx2y1z2 = xix1*yiy2*ziz1 * x2x1*y1y2*z2z1 * valin(kx+1,ky,kz+1)
    facx2y2z1 = xix1*yiy1*ziz2 * x2x1*y2y1*z1z2 * valin(kx+1,ky+1,kz)
    facx2y2z2 = xix1*yiy1*ziz1 * x2x1*y2y1*z2z1 * valin(kx+1,ky+1,kz+1)

    !*** Final value
    valout = facx1y1z1 + facx1y1z2 + facx1y2z1 + facx1y2z2 + facx2y1z1 + facx2y1z2 + facx2y2z1 + facx2y2z2
  end function interp3_0_dp

  pure subroutine append_i(list, value)
    integer(kind = IPRE), dimension(:), allocatable, intent(inout) :: list
    integer(kind = IPRE), intent(in) :: value
    integer(kind = IPRE), dimension(:), allocatable :: tmp

    if (.not. allocated(list)) then
      allocate(list(1))
      list(1) = value
      return
    end if
    allocate(tmp(size(list) + 1))
    tmp = [ list, value ]
    call move_alloc(tmp, list)
    return
  end subroutine append_i

  pure subroutine append_r8(list, value)
    real(kind = DPRE), dimension(:), allocatable, intent(inout) :: list
    real(kind = DPRE), intent(in) :: value
    real(kind = DPRE), dimension(:), allocatable :: tmp

    if (.not. allocated(list)) then
      allocate(list(1))
      list(1) = value
      return
    end if
    allocate(tmp(size(list) + 1))
    tmp = [ list, value ]
    call move_alloc(tmp, list)
    return
  end subroutine append_r8

  pure subroutine append_ch_name(list, value)
    character(len = MAX_NAME_LEN), dimension(:), allocatable, intent(inout) :: list
    character(len = MAX_NAME_LEN), intent(in) :: value
    character(len = MAX_NAME_LEN), dimension(:), allocatable :: tmp

    if (.not. allocated(list)) then
      allocate(list(1))
      list(1) = value
      return
    end if
    allocate(tmp(size(list) + 1))
    tmp = [ list, value ]
    call move_alloc(tmp, list)
    return
  end subroutine append_ch_name

  pure function find_str_loc(strings, value) result(loc)
    character(len=*), dimension(:), intent(in) :: strings
    character(len=*), intent(in) :: value
    integer :: loc, i
    
    loc = 0
    do i = 1, size(strings)
      if (trim(strings(i)) == trim(value)) then
        loc = i
        exit
      end if
    end do
  end function find_str_loc

  pure function find_int_loc(ints, value) result(loc)
    integer, dimension(:), intent(in) :: ints
    integer, intent(in) :: value
    integer :: loc, i
    
    loc = 0
    do i = 1, size(ints)
      if (ints(i) == value) then
        loc = i
        exit
      end if
    end do
  end function find_int_loc

  subroutine gradient_2_geo(f, lon, lat, tx, ty)
    ! argument type, intent(inout) :: gradient_2_geo
    real(kind=dp), dimension(:,:), intent(in) :: f
    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: tx, ty
    real(kind=dp), dimension(:,:), allocatable :: dx, dy
    real(kind=dp), dimension(:), allocatable :: dx1, dy1
    integer :: nx, ny, ix, iy
    real(kind=dp), dimension(:,:), allocatable :: gfx, gbx, gfy, gby

    nx = size(f,1)
    ny = size(f,2)
    dx1 = zeros(nx)
    dy1 = zeros(ny)
    tx = zeros(nx, ny)
    ty = zeros(nx, ny)
    gfx = zeros(nx, ny)
    gbx = zeros(nx, ny)
    gfy = zeros(nx, ny)
    gby = zeros(nx, ny)

    dx1(2:) = diff(lon)
    dy1(2:) = diff(lat)
    call meshgrid_ij(dx1, dy1, dx, dy)
    dx = radius*dx*deg2rad
    dy = radius*dy*deg2rad
    
    do iy=1,ny
      dx(:,iy) = dx(:, iy) * cosd(lat(iy))
      do ix = 1, nx
        if (iy>1) gfy(ix, iy) = f(ix,iy)-f(ix,iy-1)
        if (iy<ny) gby(ix,iy) = f(ix,iy+1)-f(ix,iy)
        if (ix>1) gfx(ix,iy) = f(ix,iy)-f(ix-1,iy)
        if (ix<nx) gbx(ix,iy) = f(ix+1,iy)-f(ix,iy)
      enddo
    enddo
    tx(2:nx-1,:) = (gfx(2:nx-1,:)+gbx(2:nx-1,:))/(2*dx(2:nx-1,:))
    tx(1,:) = gbx(1,:)/dx(2, :)
    tx(nx,:) = gfx(nx,:)/dx(nx,:)
    ty(:,2:ny-1) = (gfy(:,2:ny-1)+gby(:,2:ny-1))/(2*dy(:,2:ny-1))
    ty(:,1) = gby(:,1)/dy(:,2)
    ty(:,ny) = gfy(:,ny)/dy(:,ny)
  end subroutine gradient_2_geo

  subroutine gradient_2(f, dx, dy, tx, ty)
    real(kind=dp), dimension(:,:), allocatable, intent(in) :: f
    real(kind=dp), intent(in) :: dx, dy
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: tx, ty
    real(kind=dp), dimension(:,:), allocatable :: gfx, gbx, gfy,gby
    integer :: nx, ny, ix, iy

    nx = size(f,1)
    ny = size(f,2)
    tx = zeros(nx, ny)
    ty = zeros(nx, ny)
    gfx = zeros(nx, ny)
    gbx = zeros(nx, ny)
    gfy = zeros(nx, ny)
    gby = zeros(nx, ny)
    do ix=1,nx
      do iy = 1, ny
        if (iy>1) gfy(ix, iy) = f(ix,iy)-f(ix,iy-1)
        if (iy<ny) gby(ix,iy) = f(ix,iy+1)-f(ix,iy)
        if (ix>1) gfx(ix,iy) = f(ix,iy)-f(ix-1,iy)
        if (ix<nx) gbx(ix,iy) = f(ix+1,iy)-f(ix,iy)
      enddo
    enddo
    tx(2:nx-1,:) = (gfx(2:nx-1,:)+gbx(2:nx-1,:))/(2*dx)
    tx(1,:) = gbx(1,:)/dx
    tx(nx,:) = gfx(nx,:)/dx
    ty(:,2:ny-1) = (gfy(:,2:ny-1)+gby(:,2:ny-1))/(2*dy)
    ty(:,1) = gby(:,1)/dy
    ty(:,ny) = gfy(:,ny)/dy   
    
  end subroutine gradient_2

  pure function gaussian_smooth_geo_2(data, lons, lats, sigma) result(smdata)
    real(kind=dp), intent(in) :: sigma ! in deg
    real(kind=dp), dimension(:), intent(in) :: lons, lats
    real(kind=dp), dimension(:,:), intent(in) :: data
    real(kind=dp), dimension(:,:), allocatable :: smdata
    real(kind=dp), dimension(:,:), allocatable :: delta, w, lala, lolo
    integer :: i, j, n1, m1, n2, m2, dims(2)
    real(kind=dp) :: sigma3, wsum, dx, dy, sigma_sq
    integer :: nx, ny

    sigma3 = sigma*3
    sigma_sq = sigma**2
    dims = shape(data)
    smdata = zeros(dims(1), dims(2))
    dx = lons(2)-lons(1)
    dy = lats(2)-lats(1)
    nx = sigma3/dx
    do j = 1, dims(2)
      ny = sigma3/cosd(lats(j))/dy
      do i = 1, dims(1)
        ! wsum = 0.
        m1 = max(1, i - nx); m2 = min(dims(1), i + nx)
        n1 = max(1, j - ny); n2 = min(dims(2), j + ny)
        call meshgrid_ij(lons(m1:m2), lats(n1:n2), lolo, lala)
        delta = gps2dist(lats(j), lons(i), lala, lolo)*km2deg
        w = exp(-delta**2 / (2.0_dp * sigma_sq))
        smdata(i, j) = smdata(i, j) + sum(w*data(m1:m2, n1:n2))
        smdata(i, j) = smdata(i, j) / sum(w)
      end do
    end do
    
  end function

  pure function smooth_1(data, z, sigma) result(smdata)
    real(kind=dp), intent(in) :: sigma
    real(kind=dp), dimension(:), intent(in) :: data, z
    real(kind=dp), dimension(:), allocatable :: smdata
    real(kind=dp), dimension(:), allocatable :: w
    integer :: i, n1, n2
    real(kind=dp) :: sigma3, wsum, sigma_sq, dz
    integer :: nx

    sigma3 = sigma*3
    sigma_sq = sigma**2
    smdata = zeros(size(data))
    dz = z(2)-z(1)
    nx = sigma3/dz
    do i = 1, size(data)
      wsum = 0.
      n1 = max(1, i - nx); n2 = min(size(data), i + nx)
      w = exp(-((arange(n1,n2) - i)*dz)**2 / (2.0_dp * sigma_sq))
      smdata(i) = sum(w*data(n1:n2))
      smdata(i) = smdata(i) / sum(w)
    end do
  end function

  pure function gps2dist_scalar(lat0, lon0, lat1, lon1) result(dist)
    ! Input parameters
    real(kind=dp), intent(in) :: lon0, lon1, lat0, lat1
    ! Local variables
    real(kind=dp) :: dlat, dlon, a, deg,dist

    dlat = (lat0 - lat1) * deg2rad
    dlon = (lon1 - lon0) * deg2rad
    a = sin(dlat * 0.5) ** 2 + sin(dlon * 0.5) ** 2 &
        * cos(lat0 * deg2rad) * cos(lat1 * deg2rad)
    deg = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    dist = radius * deg
  end function gps2dist_scalar

  pure function gps2dist_1(lat0, lon0, lat1, lon1) result(dist)
    ! Input parameters
    real(kind=dp), intent(in) :: lon0, lat0
    real(kind=dp), dimension(:), intent(in) :: lon1, lat1
    real(kind=dp), dimension(:), allocatable :: dist
    ! Local variables
    real(kind=dp), dimension(:), allocatable :: dlat, dlon, a, deg

    dlat = (lat0 - lat1) * deg2rad
    dlon = (lon1 - lon0) * deg2rad
    a = sin(dlat * 0.5) ** 2 + sin(dlon * 0.5) ** 2 &
        * cos(lat0 * deg2rad) * cos(lat1 * deg2rad)
    deg = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    dist = radius * deg
  end function gps2dist_1

  pure function gps2dist_2(lat0, lon0, lat1, lon1) result(dist)
    ! Input parameters
    real(kind=dp), intent(in) :: lon0, lat0
    real(kind=dp), dimension(:,:), intent(in) :: lon1, lat1
    real(kind=dp), dimension(:,:), allocatable :: dist
    ! Local variables
    real(kind=dp), dimension(:,:), allocatable :: dlat, dlon, a, deg

    dlat = (lat0 - lat1) * deg2rad
    dlon = (lon1 - lon0) * deg2rad
    a = sin(dlat * 0.5) ** 2 + sin(dlon * 0.5) ** 2 &
        * cos(lat0 * deg2rad) * cos(lat1 * deg2rad)
    deg = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    dist = radius * deg
  end function gps2dist_2

! Locate with bisection
  subroutine locate_bissection(valx,n,x,ind)

    integer, intent(in) :: n
    real(kind=DPRE), intent(in) :: x

    real(kind=DPRE), dimension(n), intent(in) :: valx

    integer, intent(out) :: ind

    integer :: jl, ju, jm

    jl = 0
    ju = n+1

    do while ( (ju-jl) > 1 )
        jm = (ju + jl) / 2
        if ( x >= valx(jm) ) then
            jl = jm
        else
            ju = jm
        endif
    enddo
    if ( x == valx(1) ) then
        ind = 1
        else if ( x == valx(n) ) then
            ind = n-1
        else
            ind = jl
    endif

  end subroutine locate_bissection
end module