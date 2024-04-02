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

program rotate_src_rec
  use shared_par
  use argparse, only: argparse_rotate_src_rec
  use src_rec_parser
  use sph2loc
  use para, ap => att_para_global
  use src_rec, sr_gr => src_rec_global_gr, sr_ph => src_rec_global_ph
  implicit none

  character(len=MAX_STRING_LEN) :: fname, outfname
  real(kind=dp) :: angle, center(2)
  type(SrcRecRaw) :: sr
  real(kind=dp), dimension(:), allocatable :: new_lat, new_lon

  call argparse_rotate_src_rec(fname, angle, center, outfname)
  call sr%read_raw_src_rec_file(fname)
  
  ! Get the center of the sources
  ! Rotate the source
  call rtp_rotation(sr%evla, sr%evlo, center(1), center(2), angle, new_lat, new_lon)
  sr%evla = new_lat; sr%evlo = new_lon
  deallocate(new_lat, new_lon)

  ! Get the center of the receivers
  ! Rotate the receiver
  call rtp_rotation(sr%stla, sr%stlo, center(1), center(2), angle, new_lat, new_lon)
  sr%stla = new_lat; sr%stlo = new_lon
  deallocate(new_lat, new_lon)

  ! Write the rotated source-receiver file
  call sr%write_raw_src_rec_file(outfname)
end program rotate_src_rec
