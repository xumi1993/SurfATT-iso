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

program surfatt_tomo
  use my_mpi
  use para, ap => att_para_global
  use src_rec, sr_gr => src_rec_global_gr, sr_ph => src_rec_global_ph
  use model, am => att_model_global
  use grid, ag_gr => att_grid_global_gr, ag_ph => att_grid_global_ph
  use measadj
  use tomo
  use setup_att_log, only: setuplog
  use argparse, only: argparse_tomo

  implicit none
  
  type(att_tomo) :: att
  character(len=MAX_STRING_LEN) :: fname
  logical :: isfwd
  real(kind=dp) :: t0, t1

  ! initialize MPI
  call init_mpi()
  call world_rank(myrank)
  call world_size(mysize)

  call cpu_time(t0)
  ! read command line arguments
  call argparse_tomo(fname, isfwd)

  ! read parameter file
  call ap%read(fname)

  ! intialize logger
  call setuplog(ap%output%log_level)

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
  call att%init(is_fwd=isfwd)

  if (isfwd) then
    ! do forward
    call att%do_forward()
  else
    ! do inversion
    call att%do_inversion()
  endif

  ! calculate CPU time 
  call cpu_time(t1)
  write(att%message, '(a,f0.2,a)') 'Elapsed CPU time: ', t1-t0, ' s'
  call write_log(att%message,1,att%module)

  ! MPI finish
  call finalize_mpi()
end program surfatt_tomo
