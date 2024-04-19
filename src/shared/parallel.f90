!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! Dimitri Komatitsch, July 2014, CNRS Marseille, France:
! added the ability to run several calculations (several earthquakes)
! in an embarrassingly-parallel fashion from within the same run;
! this can be useful when using a very large supercomputer to compute
! many earthquakes in a catalog, in which case it can be better from
! a batch job submission point of view to start fewer and much larger jobs,
! each of them computing several earthquakes in parallel.
! To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1 in the Par_file.
! To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,
! each of them being labeled "my_local_mpi_comm_world", and we use them
! in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case
! we need to kill the entire run.
! When that option is on, of course the number of processor cores used to start
! the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,
! all the individual runs must use the same number of processor cores,
! which as usual is NPROC in the input file DATA/Par_file,
! and thus the total number of processor cores to request from the batch system
! should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.
! All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on
! (with exactly four digits).

!-------------------------------------------------------------------------------------------------
!
! Parallel routines.  All MPI calls belong in this file!
!
!-------------------------------------------------------------------------------------------------

module my_mpi

! main parameter module for specfem simulations

  use shared_par
  use mpi

  implicit none

  integer :: my_local_mpi_comm_world, my_node_mpi_comm_world, my_local_mpi_comm_for_bcast

  interface bcast_all
    module procedure bcast_all_ch_array, bcast_all_cr, bcast_all_dp, bcast_all_i,&
                     bcast_all_l_array, bcast_all_singlecr, bcast_all_dp_2, bcast_all_dp_3,&
                     bcast_all_singledp, bcast_all_singlei, bcast_all_singlel, &
                     bcast_all_string
  end interface

  interface min_all
    module procedure min_all_cr, min_all_dp, min_all_i
  end interface

  interface max_all
    module procedure max_all_cr, max_all_dp, max_all_i
  end interface

  interface max_all_all
    module procedure max_all_all_cr, max_all_all_dp, max_all_all_i
  end interface

  interface sum_all
    module procedure sum_all_1Darray_dp, sum_all_2Darray_dp, &
                     sum_all_3Darray_dp, sum_all_4Darray_dp, &
                     sum_all_1Darray_i, sum_all_2Darray_i,sum_all_3Darray_i,&
                     sum_all_cr, sum_all_dp, sum_all_i
  end interface

  interface sum_all_all
    module procedure sum_all_all_1Darray_dp, sum_all_all_2Darray_dp,&
                     sum_all_all_3Darray_dp, sum_all_all_4Darray_dp,&
                     sum_all_all_cr, sum_all_all_i
  end interface

  interface send
   module procedure send_dp, send_i, send_ch_array
  end interface

  interface recv
    module procedure recv_dp, recv_i, recv_ch_array
  end interface

  interface sync_from_main_rank
    module procedure sync_from_main_rank_ch, sync_from_main_rank_dp, sync_from_main_rank_i,&
                     sync_from_main_rank_dp_2, sync_from_main_rank_dp_3, sync_from_main_rank_dp_4
  end interface

  contains
!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  integer :: ier
  integer, dimension(:,:), allocatable :: rank_map_loc

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read before calling world_split()
  ! thus read the parameter file
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  ! if (myrank == 0) then
  !   call open_parameter_file_from_master_only(ier)
  !   ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read
  !   call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
  !   if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
  !   call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
  !   if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'
  !   ! close parameter file
  !   call close_parameter_file()
  ! endif

  ! broadcast parameters read from master to all processes
  my_local_mpi_comm_world = MPI_COMM_WORLD
  ! call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
  ! call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)
  ! my_local_mpi_comm_for_bcast = MPI_COMM_NULL

! create sub-communicators if needed, if running more than one earthquake from the same job
  ! call world_split()
  call world_rank(myrank)
  call world_size(mysize)

  call MPI_Comm_split_type(my_local_mpi_comm_world, MPI_COMM_TYPE_SHARED, myrank, &
                           MPI_INFO_NULL, my_node_mpi_comm_world, ier)
  call MPI_Comm_rank(my_node_mpi_comm_world, local_rank, ier)
  call MPI_Comm_size(my_node_mpi_comm_world, local_size, ier)

  allocate(rank_map_loc(mysize, 2), rank_map(mysize, 2))
  rank_map_loc = 0

  rank_map_loc(myrank+1, 1) = myrank
  rank_map_loc(myrank+1, 2) = local_rank
  call sum_all(rank_map_loc, rank_map, mysize, 2)
  
  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  ! use my_mpi

  implicit none

  integer :: ier


  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
! !

  subroutine abort_mpi()

  ! use my_mpi
  use constants, only: MAX_STRING_LEN
  ! use shared_input_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: my_local_rank,my_global_rank,ier
  logical :: run_file_exists
  character(len=MAX_STRING_LEN) :: filename

  ! get my local rank and my global rank (in the case of simultaneous jobs, for which we split
  ! the MPI communicator, they will be different; otherwise they are the same)
  call world_rank(my_local_rank)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_global_rank,ier)

  ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'error, program ended in exit_MPI'

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

  ! use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_world,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_all

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine synchronize_all_comm(comm)
!  end subroutine synchronize_all_comm

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

  ! use my_mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

!  integer function null_process()
!  end function null_process

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine test_request(request,flag_result_test)
!  end subroutine test_request

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_req(req)

  ! use my_mpi

  implicit none

  integer :: req

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

  end subroutine wait_req


!-------------------------------------------------------------------------------------------------
!
! MPI broadcasting helper
!
!-------------------------------------------------------------------------------------------------

  subroutine bcast_all_i(buffer, countval)

  ! use my_mpi

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  ! use my_mpi

  implicit none

  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlel(buffer)

  ! use my_mpi

  implicit none

  logical :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr(buffer, countval)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: countval
  real(kind=cr), dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlecr(buffer)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlecr

!
!-------------------------------------------------------------------------------------------------
!

  ! subroutine bcast_all_r(buffer, countval)

  ! ! use my_mpi

  ! implicit none

  ! integer :: countval
  ! real, dimension(countval) :: buffer

  ! integer :: ier

  ! call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

  ! end subroutine bcast_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp(buffer, countval)

  ! use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp

  subroutine bcast_all_dp_2(buffer, countval)

  ! use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(:,:) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp_2

subroutine bcast_all_dp_3(buffer, countval)

  ! use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(:,:,:) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp_3
!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

  ! use my_mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array(buffer,countval,STRING_LEN)

    ! use my_mpi

    implicit none

    integer :: countval, STRING_LEN

    character(len=STRING_LEN), dimension(countval) :: buffer

    integer :: ier

    call MPI_BCAST(buffer,STRING_LEN*countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_ch_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l_array(buffer, countval)

    ! use my_mpi
    implicit none
    integer    :: countval
    logical, dimension(countval) :: buffer
    integer :: ier

    call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_l_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_string(buffer)

  ! use my_mpi
  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_string

!
!---- broadcast to send the mesh and model to other simultaneous runs
!

  ! subroutine bcast_all_i_for_database(buffer, countval)

  ! ! use my_mpi
  ! ! use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  ! implicit none

  ! integer countval
  ! ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  ! integer :: buffer

  ! integer ier

  ! if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  ! call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  ! end subroutine bcast_all_i_for_database

!
!-------------------------------------------------------------------------------------------------
!

  ! subroutine bcast_all_l_for_database(buffer, countval)

  ! ! use my_mpi
  ! use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  ! implicit none

  ! integer countval
  ! ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  ! logical :: buffer

  ! integer ier

  ! if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  ! call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  ! end subroutine bcast_all_l_for_database

!
!-------------------------------------------------------------------------------------------------
!

!   subroutine bcast_all_cr_for_database(buffer, countval)

!   ! use my_mpi
!   use constants, only: cr
!   use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

!   implicit none

  

!   integer countval
!   ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
!   ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
!   ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
!   real(kind=cr) :: buffer

!   integer ier

!   if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

!   call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)

!   end subroutine bcast_all_cr_for_database

! !
! !-------------------------------------------------------------------------------------------------
! !

!   subroutine bcast_all_dp_for_database(buffer, countval)

!   ! use my_mpi
!   use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

!   implicit none

!   integer countval
!   ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
!   ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
!   ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
!   double precision :: buffer

!   integer ier

!   if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

!   call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_for_bcast,ier)

!   end subroutine bcast_all_dp_for_database

! !
! !-------------------------------------------------------------------------------------------------
! !

!   subroutine bcast_all_r_for_database(buffer, countval)

!   ! use my_mpi
!   use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

!   implicit none

!   integer countval
!   ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
!   ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
!   ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
!   real :: buffer

!   integer ier

!   if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

!   call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)

!   end subroutine bcast_all_r_for_database

!-------------------------------------------------------------------------------------------------
!
! MPI math helper
!
!-------------------------------------------------------------------------------------------------

  subroutine min_all_i(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_i

  subroutine max_all_all_i(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_allreduce_i(buffer,countval)

  ! use my_mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)
  if (ier /= 0) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr):: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_MIN,my_local_mpi_comm_world,ier)

  end subroutine min_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine max_allreduce_cr(sendbuf, recvbuf)
!  end subroutine max_allreduce_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision, dimension(2) :: sendbuf,recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,my_local_mpi_comm_world,ier)

  end subroutine maxloc_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine any_all_l(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  logical :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,my_local_mpi_comm_world,ier)

  end subroutine any_all_l

!
!-------------------------------------------------------------------------------------------------
!

! MPI summations

  subroutine sum_all_i(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_i(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  real(kind=cr) sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_REAL,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_dp


  subroutine sum_all_all_dp(sendbuf, recvbuf)

  ! use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)


  end subroutine sum_all_all_dp
!
!-------------------------------------------------------------------------------------------------
!

 subroutine sum_all_1Darray_i(sendbuf, recvbuf, nx)

  ! use my_mpi

  implicit none

  integer :: nx
  integer, dimension(nx) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_i

   subroutine sum_all_2Darray_i(sendbuf, recvbuf, nx, ny)

  ! use my_mpi

  implicit none

  integer :: nx, ny
  integer, dimension(nx,ny) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx*ny,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_2Darray_i

  subroutine sum_all_3Darray_i(sendbuf, recvbuf, nx, ny, nz)

  ! use my_mpi

  implicit none

  integer :: nx, ny, nz
  integer, dimension(nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_3Darray_i

  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)

  ! use my_mpi

  implicit none

  integer :: nx
  double precision, dimension(nx) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_dp

  subroutine sum_all_all_1Darray_dp(sendbuf, recvbuf, nx)

    implicit none

    integer :: nx
    double precision, dimension(nx) :: sendbuf, recvbuf
    integer :: ier

    call MPI_ALLREDUCE(sendbuf,recvbuf,nx,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_1Darray_dp
!
!-------------------------------------------------------------------------------------------------
!
 subroutine sum_all_all_2Darray_dp(sendbuf, recvbuf, nx,ny)
  integer :: nx,ny
  double precision, dimension(nx,ny) :: sendbuf, recvbuf
  integer :: ier

  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_all_2Darray_dp

 subroutine sum_all_2Darray_dp(sendbuf, recvbuf, nx,ny)
  integer :: nx,ny
  double precision, dimension(nx,ny) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_2Darray_dp


 subroutine sum_all_all_3Darray_dp(sendbuf, recvbuf, nx,ny,nz)
  integer :: nx,ny,nz
  double precision, dimension(nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_all_3Darray_dp

 subroutine sum_all_3Darray_dp(sendbuf, recvbuf, nx,ny,nz)
  integer :: nx,ny,nz
  double precision, dimension(nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_3Darray_dp

 subroutine sum_all_4Darray_dp(sendbuf, recvbuf, nr,nx,ny,nz)
  integer :: nx,ny,nz, nr
  double precision, dimension(nr,nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nr*nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nr*nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_4Darray_dp


 subroutine sum_all_all_4Darray_dp(sendbuf, recvbuf, nr,nx,ny,nz)
  integer :: nx,ny,nz,nr
  double precision, dimension(nr,nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  ! call MPI_ALLREDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)
  call MPI_ALLREDUCE(sendbuf,recvbuf,nr*nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

 end subroutine sum_all_all_4Darray_dp


!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------

! asynchronuous send/receive

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: sendcount, dest, sendtag, req
  real(kind=cr), dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,MPI_REAL,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)

  ! use my_mpi

  implicit none

  integer :: sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_i

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine isend_dp(sendbuf, sendcount, dest, sendtag, req)
!  end subroutine isend_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: recvcount, dest, recvtag, req
  real(kind=cr), dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,MPI_REAL,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine irecv_dp(recvbuf, recvcount, dest, recvtag, req)
!  end subroutine irecv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  ! use my_mpi

  implicit none

  integer :: recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf
  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_i

!
!-------------------------------------------------------------------------------------------------
!

! synchronuous/blocking send/receive

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  ! use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  ! use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: recvcount,dest,recvtag
  real(kind=cr),dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_REAL,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recvv_cr

  subroutine recv_ch_array(recvbuf, recvcount, nlen, dest, recvtag)

  ! use my_mpi

  implicit none

  integer :: dest,recvtag,nlen
  integer :: recvcount
  character(len=nlen),dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_CHARACTER,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_ch_array
!
!-------------------------------------------------------------------------------------------------
!

!  subroutine recv_singlei(recvbuf, dest, recvtag)
!  end subroutine recv_singlei

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine recv_singlel(recvbuf, dest, recvtag)
!  end subroutine recv_singlel

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine send_ch(sendbuf, sendcount, dest, sendtag)
!  end subroutine send_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  ! use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  ! use my_mpi

  implicit none

  integer :: dest,sendcount,ier
  integer :: tag = 100
  integer, dimension(sendcount) :: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,tag,my_local_mpi_comm_world,ier)

  end subroutine send_i_t

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_i_t(recvbuf,recvcount,source)

  ! use my_mpi

  implicit none

  integer :: source,recvcount,ier
  integer :: tag = 100
  integer, dimension(recvcount) :: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,source,tag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i_t

!
!-------------------------------------------------------------------------------------------------
!
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  ! use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_dp


  subroutine send_ch_array(sendbuf, sendcount, nlen, dest, sendtag)

  ! use my_mpi

  implicit none

  integer :: dest,sendtag,nlen
  integer :: sendcount
  character(len=nlen),dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_CHARACTER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_ch_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer sendcount,dest,sendtag
  real(kind=cr),dimension(sendcount) :: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_REAL,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine sendv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sendrecv_cr(sendbuf, sendcount, dest, sendtag, &
!                         recvbuf, recvcount, source, recvtag)
!  end subroutine sendrecv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sendrecv_dp(sendbuf, sendcount, dest, sendtag, &
!                         recvbuf, recvcount, source, recvtag)
!  end subroutine sendrecv_dp


!-------------------------------------------------------------------------------------------------
!
! MPI gather helper
!
!-------------------------------------------------------------------------------------------------

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  ! use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  ! use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,sendcnt,MPI_INTEGER, &
                     recvbuf,recvcount,MPI_INTEGER, &
                     my_local_mpi_comm_world,ier)

  end subroutine gather_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  ! use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: sendcnt, recvcount, NPROC
  real(kind=cr), dimension(sendcnt) :: sendbuf
  real(kind=cr), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_REAL, &
                  recvbuf,recvcount,MPI_REAL, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  ! use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_dp

!
!-------------------------------------------------------------------------------------------------
!

! unused so far...
!
!  subroutine gather_all_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)
!
!  ! use my_mpi
!
!  implicit none
!
!  integer :: sendcnt, recvcount, NPROC
!  double precision, dimension(sendcnt) :: sendbuf
!  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf
!
!  integer :: ier
!
!  call MPI_ALLGATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
!                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
!                  my_local_mpi_comm_world,ier)
!
!  end subroutine gather_all_all_dp
!

!
!-------------------------------------------------------------------------------------------------
!


  subroutine gather_all_all_cr(sendbuf, recvbuf, counts, NPROC)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer NPROC,counts
  real(kind=cr), dimension(counts) :: sendbuf
  real(kind=cr), dimension(counts,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_ALLGATHER(sendbuf,counts,MPI_REAL,recvbuf,counts,MPI_REAL,my_local_mpi_comm_world,ier)

  end subroutine gather_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  ! use my_mpi

  implicit none

  

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_INTEGER,recvbuf,recvcount,recvoffset,MPI_INTEGER, &
                   0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_i


!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  ! use my_mpi
  use constants, only: cr

  implicit none

  

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=cr), dimension(sendcnt) :: sendbuf
  real(kind=cr), dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_REAL, &
                  recvbuf,recvcount,recvoffset,MPI_REAL, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gatherv_all_ch(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, dim2, NPROC)

  use constants
  ! use my_mpi

  implicit none

  integer :: sendcnt, dim1, dim2, NPROC

  character(len=dim2), dimension(NPROC) :: sendbuf
  character(len=dim2), dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer :: ier

  call MPI_Allgatherv(sendbuf,sendcnt,MPI_CHARACTER, &
                  recvbuf,recvcnt,recvoffset,MPI_CHARACTER, &
                  my_local_mpi_comm_world,ier)

  end subroutine all_gatherv_all_ch


! not used yet...
!  subroutine gatherv_all_dp(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)
!
!  ! use my_mpi
!  use constants, only: cr
!
!  implicit none
!
!  
!
!  integer :: sendcnt,recvcounttot,NPROC
!  integer, dimension(NPROC) :: recvcount,recvoffset
!  double precision, dimension(sendcnt) :: sendbuf
!  double precision, dimension(recvcounttot) :: recvbuf
!
!  integer :: ier
!
!  call MPI_GATHERV(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
!                  recvbuf,recvcount,recvoffset,MPI_DOUBLE_PRECISION, &
!                  0,my_local_mpi_comm_world,ier)
!
!  end subroutine gatherv_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine gatherv_all_r(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)
!  end subroutine gatherv_all_r

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine scatter_all_singlei(sendbuf, recvbuf, NPROC)
!  end subroutine scatter_all_singlei


!-------------------------------------------------------------------------------------------------
!
! MPI world helper
!
!-------------------------------------------------------------------------------------------------

  subroutine world_size(sizeval)

  ! use my_mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  ! use my_mpi

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_duplicate(comm)

  ! use my_mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  call MPI_COMM_DUP(my_local_mpi_comm_world,comm,ier)
  if (ier /= 0) stop 'error duplicating my_local_mpi_comm_world communicator'

  end subroutine world_duplicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm(comm)

  ! use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = my_local_mpi_comm_world

  end subroutine world_get_comm

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine world_get_comm_self(comm)
!  end subroutine world_get_comm_self

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine world_get_info_null(info)
!  end subroutine world_get_info_null


  subroutine exit_MPI(myrank,error_msg)

    use constants

    implicit none

  ! identifier for error message file
    integer, parameter :: IERROR = 30

    integer myrank
    character(len=*) :: error_msg

    character(len=MAX_STRING_LEN) :: outputname

  ! write error message to screen
    write(*,*) error_msg(1:len(error_msg))
    write(*,*) 'Error detected, aborting MPI... proc ',myrank
    call abort_mpi()
  end subroutine

  subroutine scatter_all_i(countval,size, rank, istart, iend)
    integer :: countval,size, rank, count, remainder, istart, iend
    ! integer, dimension(:), allocatable, intent(out) :: local_index

    if (countval < size) then  
        if (rank < countval) then
            istart = rank + 1
            iend = istart
        else
            istart = 1  ! invalid index
            iend = 0    ! ensures empty range
        end if
    else
        count = countval / size        ! Number of elements per rank
        remainder = mod(countval, size)           ! If n is not a perfect multiple of size

        if (rank < remainder) then
            istart = rank * (count + 1) + 1
            iend = istart + count
        else
            istart = rank * count + remainder + 1
            iend = istart + count - 1 
        end if
    end if


  end subroutine scatter_all_i

  ! 
  subroutine sync_from_main_rank_dp(buffer, countval)
    integer, intent(in) :: countval
    double precision, dimension(:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send(buffer, countval, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call recv(buffer, countval, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp


  subroutine sync_from_main_rank_dp_2(buffer, nx, ny)
    integer, intent(in) :: nx, ny
    double precision, dimension(:,:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send_dp(buffer, nx*ny, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call recv_dp(buffer, nx*ny, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp_2


  subroutine sync_from_main_rank_dp_3(buffer, nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(:,:,:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send_dp(buffer, nx*ny*nz, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call recv_dp(buffer, nx*ny*nz, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp_3


  subroutine sync_from_main_rank_dp_4(buffer, n1, n2, n3, n4)
    integer, intent(in) :: n1, n2, n3, n4
    double precision, dimension(:,:,:,:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send_dp(buffer, n1*n2*n3*n4, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call recv_dp(buffer, n1*n2*n3*n4, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp_4


  subroutine sync_from_main_rank_i(buffer, countval)
    integer, intent(in) :: countval
    integer, dimension(:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send_i(buffer, countval, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call send_i(buffer, countval, 0, tag)
    endif
  end subroutine sync_from_main_rank_i

  subroutine sync_from_main_rank_ch(buffer, countval, nlen)
    integer, intent(in) :: countval, nlen
    character(len=nlen), dimension(:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (myrank == 0) then
      do i = 2, mysize
        if (rank_map(i, 2) == 0) then
          call send_ch_array(buffer, countval, nlen, rank_map(i, 1), tag)
        endif
      enddo
    elseif (local_rank == 0) then
      call recv_ch_array(buffer, countval, nlen, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_ch

  ! -------------------------------------------------------------------------------------------------
  subroutine prepare_shm_array_dp_1d(buffer, n_elem, win)
    USE, INTRINSIC :: ISO_C_BINDING
    double precision, dimension(:), pointer :: buffer
    integer :: ierr, istat,n_elem,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    ! call world_rank(myrank)
    n = n_elem
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE = [n_elem])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_1d

  subroutine prepare_shm_array_dp_2d(buffer, nx, ny, win)
    USE, INTRINSIC :: ISO_C_BINDING
    double precision, dimension(:,:), pointer :: buffer
    integer :: ierr,nx,ny,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    ! call world_rank(myrank)
    n = nx*ny
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE=[nx, ny])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_2d

  subroutine prepare_shm_array_i_2d(buffer, nx, ny, win)
    USE, INTRINSIC :: ISO_C_BINDING
    integer, dimension(:,:), pointer :: buffer
    integer :: ierr,nx,ny,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, int_size
    type(C_PTR) :: c_window_ptr
    
    n = nx*ny
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_INTEGER, int_size, ierr)
    size = n * int_size
    call MPI_Win_allocate_shared(size, int_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, int_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE=[nx, ny])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_i_2d

  subroutine prepare_shm_array_dp_3d(buffer, nx, ny, nz, win)
    USE, INTRINSIC :: ISO_C_BINDING
    double precision, dimension(:,:,:), pointer :: buffer
    integer :: ierr,nx,ny,nz,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    n = nx*ny*nz
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE=[nx, ny, nz])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_3d

  subroutine prepare_shm_array_dp_4d(buffer, n1, n2, n3, n4, win)
    USE, INTRINSIC :: ISO_C_BINDING
    double precision, dimension(:,:,:,:), pointer :: buffer
    integer :: ierr,n1,n2,n3,n4,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    n = n1*n2*n3*n4
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE=[n1,n2,n3,n4])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_4d

  subroutine prepare_shm_array_ch_1d(buffer, n_elem, nlen, win)
    USE, INTRINSIC :: ISO_C_BINDING
    integer :: ierr, istat,n_elem,n,nlen
    character(len=nlen), dimension(:), pointer :: buffer
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, char_size
    type(C_PTR) :: c_window_ptr
    
    n = n_elem*nlen
    if(local_rank /= 0) n = 0
    CALL MPI_Type_size(MPI_CHARACTER, char_size, ierr)
    size = n * char_size
    call MPI_Win_allocate_shared(size, char_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (local_rank /= 0) then
      call MPI_Win_shared_query(win, 0, size, char_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE = [n_elem])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_ch_1d

end module my_mpi
