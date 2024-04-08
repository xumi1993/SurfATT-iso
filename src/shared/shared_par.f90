module constants
  include "constants.h"


  ! create a copy of the original output file path, to which we may add a "run0001/", "run0002/", "run0003/" prefix later
  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES = 'OUTPUT_FILES'

  ! if doing simultaneous runs for the same mesh and model, see who should read the mesh and the model and broadcast it to others
  ! we put a default value here
  real(kind=dp), parameter :: pi = 3.14159265359
  real(kind=dp), parameter :: radius = 6371.0d0
  real(kind=dp), parameter :: deg2rad = pi/180.0d0
  real(kind=dp), parameter :: rad2deg = 180.0d0/pi
  real(kind=dp), parameter :: km2deg = 1/(6371.0d0*pi/180.0d0)

  character(len=MAX_STRING_LEN),parameter :: srfile = 'src_rec_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: modfile = 'model_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: log_fname='output_attsurf_tomo.log'


  real(kind=dp), parameter :: precond_thres = 1.0d-2 
  integer, parameter :: iter_store = 6

end module constants

module shared_par
  use constants
  implicit none

  integer :: myrank, mysize
  integer :: LID
  integer :: loglevel
  real(kind=dp), dimension(:,:,:), allocatable :: gradient_s
end module shared_par