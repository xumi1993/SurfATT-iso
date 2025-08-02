module constants
  include "constants.h"
  include "version.h"

  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES = 'OUTPUT_FILES'
  integer, parameter :: GROUP_VEL = 1, PHASE_VEL = 0
  integer, parameter :: RAYLEIGH = 2, LOVE = 1
  
  ! we put a default value here
  real(kind=dp), parameter :: pi = 3.14159265359
  real(kind=dp), parameter :: radius = 6371.0d0
  real(kind=dp), parameter :: deg2rad = pi/180.0d0
  real(kind=dp), parameter :: rad2deg = 180.0d0/pi
  real(kind=dp), parameter :: km2deg = 1/(6371.0d0*pi/180.0d0)
  real(kind=dp), parameter :: rho_scaling = 0.33

  character(len=MAX_STRING_LEN),parameter :: srfile = 'src_rec_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: modfile = 'model_iter.h5'
  character(len=MAX_STRING_LEN),parameter :: log_basename='output_attsurf_tomo.log'


  real(kind=dp), parameter :: precond_thres = 1.0d-4
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
  real(kind=dp), dimension(:,:,:,:), allocatable :: gradient_s, direction
  character(len=MAX_STRING_LEN) :: log_fname, model_fname
  character(len=MAX_STRING_LEN) :: VERSION
  integer :: m_store = 5
  integer :: nker = 1

contains
  subroutine get_version()
    write(VERSION, '(I0,".",I0,".",I0)') PROJECT_VERSION_MAJOR, PROJECT_VERSION_MINOR, PROJECT_VERSION_PATCH
  end subroutine get_version
end module shared_par