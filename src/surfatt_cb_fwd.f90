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

program surfatt_cb_fwd
  use my_mpi
  use para, ap => att_para_global
  use src_rec, sr_gr => src_rec_global_gr, sr_ph => src_rec_global_ph
  use model, am => att_model_global
  use measadj
  use tomo
  use setup_att_log, only: setuplog
  use argparse, only: argparse_cb_fwd
  ! use stdlib_io_npy, only: save_npy

  implicit none
  
  type(att_tomo) :: att
  character(len=MAX_STRING_LEN) :: fname
  integer, dimension(3) :: ncb
  real(kind=dp) :: pert, hmarg, anom_size, max_noise
  logical :: only_vs

  ! initialize MPI
  call init_mpi()

  ! read command line arguments
  call argparse_cb_fwd(fname, ncb, pert, hmarg, anom_size, max_noise, only_vs)

  ! read parameter file
  call ap%read(fname)

  ! intialize logger
  call setuplog(ap%output%log_level)

  ! read dispersion data
  if (ap%data%vel_type(1)) call sr_ph%read(type=0)
  if (ap%data%vel_type(2)) call sr_gr%read(type=1)
  call merge_sta()

  ! initialize model type
  call am%init()

  ! construct initial model
  call am%get_init_model()

  ! add perturbations
  call am%add_pert(ncb(1),ncb(2),ncb(3), pert_vel=pert,&
                   hmarg=hmarg, anom_size=anom_size,only_vs=only_vs)

  ! initialize grid
  if (ap%data%vel_type(1)) then
    call ag_ph%init(1)
    call ag_ph%get_topo()
  endif
  if (ap%data%vel_type(2)) then
    call ag_gr%init(2)
    call ag_gr%get_topo()
  endif

  call att%init(is_fwd=.true.)
  call am%write('target_model')

  ! do forward
  call att%do_forward(max_noise=max_noise)

  ! MPI finish
  call finalize_mpi()
end program surfatt_cb_fwd
