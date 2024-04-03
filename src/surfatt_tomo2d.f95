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

program surfatt_tomo2d
  use my_mpi
  use para, ap => att_para_global
  use src_rec, sr_gr => src_rec_global_gr, sr_ph => src_rec_global_ph
  use model, am => att_model_global
  use grid, ag_gr => att_grid_global_gr, ag_ph => att_grid_global_ph
  use measadj
  use tomo2d
  use setup_att_log, only: setuplog
  use argparse, only: argparse_tomo2d
  ! use stdlib_io_npy, only: save_npy

  implicit none
  
  type(att_tomo_2d) :: att
  character(len=MAX_STRING_LEN) :: fname
  logical :: isfwd
  integer, dimension(2) :: ncb
  real(kind=dp) :: pert_vel, hmarg

  ! initialize MPI
  call init_mpi()
  call world_rank(myrank)
  call world_size(mysize)

  ! read command line arguments
  call argparse_tomo2d(fname, isfwd, ncb, pert_vel, hmarg)

  ! read parameter file
  call ap%read(fname)

  ! intialize logger
  call setuplog()

  ! read dispersion data
  if (ap%data%vel_type(1)) call sr_ph%read(0)
  if (ap%data%vel_type(2)) call sr_gr%read(1)
  call merge_sta()

  ! initialize model type
  call am%init()

  ! construct initial model
  call am%get_init_model()

  ! initialize grid
  if (ap%data%vel_type(1)) then
    call ag_ph%init(1)
    call ag_ph%get_topo()
  endif
  if (ap%data%vel_type(2)) then
    call ag_gr%init(2)
    call ag_gr%get_topo()
  endif

  ! initial inverison
  call att%init()

  if (isfwd) then
    ! do forward
    call att%do_forward(ncb, pert_vel, hmarg)
    call am%write('target_model')
  else
    ! do inversion
    call att%do_inversion()
  endif

  ! MPI finish
  call finalize_mpi()
end program surfatt_tomo2d
