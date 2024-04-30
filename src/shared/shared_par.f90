module constants
  include "constants.h"

  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES = 'OUTPUT_FILES'
  
  ! we put a default value here
  real(kind=dp), parameter :: pi = 3.14159265359
  real(kind=dp), parameter :: radius = 6371.0d0
  real(kind=dp), parameter :: deg2rad = pi/180.0d0
  real(kind=dp), parameter :: rad2deg = 180.0d0/pi
  real(kind=dp), parameter :: km2deg = 1/(6371.0d0*pi/180.0d0)

  character(len=MAX_STRING_LEN),parameter :: srfile = 'src_rec_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: modfile = 'model_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: log_basename='output_attsurf_tomo.log'


  real(kind=dp), parameter :: precond_thres = 1.0d-4
  integer, parameter :: max_iter_1d = 40
  integer, parameter :: m_store = 5
  integer, parameter :: iter_start = 0

end module constants

module shared_par
  use constants
  implicit none

  integer :: myrank, mysize
  integer :: local_rank, local_size
  integer, dimension(:,:), allocatable :: rank_map
  integer :: LID
  integer :: loglevel
  real(kind=dp), dimension(:,:,:), allocatable :: gradient_s, direction
  character(len=MAX_STRING_LEN) :: log_fname, model_fname
end module shared_par