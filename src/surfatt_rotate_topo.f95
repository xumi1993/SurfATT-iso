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

program rotate_topo
  use shared_par
  use argparse, only: argparse_rotate_topo
  use src_rec_parser
  use sph2loc
  use topo
  implicit none

  character(len=MAX_STRING_LEN) :: fname, outfname
  real(kind=dp) :: angle, center(2), xrange(2), yrange(2)
  real(kind=dp), dimension(:), allocatable :: new_lat, new_lon
  type(att_topo) :: at

  ! Parse command line arguments
  call argparse_rotate_topo(fname, xrange, yrange, angle, center, outfname)

  ! Read in the topography file
  call at%read(fname)

  ! Rotate the topography
  call at%rotate(xrange(1), xrange(2), yrange(1), yrange(2), center(1), center(2), angle)

  ! Write the rotated topography to a new file
  call at%write(outfname)
end program rotate_topo
