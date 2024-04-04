!> hdf5_interface
module hdf5_interface
  use H5LT
  implicit none

  public :: hdf5_file

  interface h5read
    module procedure h5read_f_1d, h5read_f_2d, h5read_f_3d
  end interface h5read
  public :: h5read

  interface h5write
    module procedure h5write_f_1d, h5write_f_2d, h5write_f_3d
  end interface h5write
  public :: h5write

  private
  
  integer, parameter :: dp = 8
  type :: hdf5_file

    character(:), allocatable :: filename
    integer(HID_T) :: lid    !< location identifier
    integer(HID_T) :: gid    !< group identifier
    integer(HID_T) :: glid   !< group location identifier
    logical :: is_init = .FALSE.


  contains
    !> open/close HDF5 file
    procedure :: open  => hdf_open_file
    procedure :: close => hdf_close_file

    !> open/close HDF5 group
    procedure :: open_group  => hdf_open_group
    procedure :: close_group => hdf_close_group

    !> HDF5 group/dataset exist status
    procedure :: exist => hdf_exist

    !> delete HDF5 group/dataset
    procedure :: delete => hdf_delete

    !> add group or dataset integer/real 0-3d
    generic   :: add => hdf_add_group,&
                        hdf_add_int,&
                        hdf_add_int1d,&
                        hdf_add_int2d,&
                        hdf_add_int3d,&
                        hdf_add_real,&
                        hdf_add_real1d,&
                        hdf_add_real2d,&
                        hdf_add_real3d,&
                        hdf_add_string

    !> get dataset integer/real 0-3d
    generic   :: get => hdf_get_int,&
                        hdf_get_int1d,&
                        hdf_get_int2d,&
                        hdf_get_int3d,&
                        hdf_get_real,&
                        hdf_get_real1d,&
                        hdf_get_real2d,&
                        hdf_get_real3d,&
                        hdf_get_string

    !> add attribute
    generic    :: adda => hdf_adda_string

    !> private methods
    procedure,private :: hdf_add_group
    procedure,private :: hdf_add_int
    procedure,private :: hdf_get_int
    procedure,private :: hdf_add_int1d
    procedure,private :: hdf_get_int1d
    procedure,private :: hdf_add_int2d
    procedure,private :: hdf_get_int2d
    procedure,private :: hdf_add_int3d
    procedure,private :: hdf_get_int3d
    procedure,private :: hdf_add_real
    procedure,private :: hdf_get_real
    procedure,private :: hdf_add_real1d
    procedure,private :: hdf_get_real1d
    procedure,private :: hdf_add_real2d
    procedure,private :: hdf_get_real2d
    procedure,private :: hdf_add_real3d
    procedure,private :: hdf_get_real3d
    procedure,private :: hdf_get_string
    procedure,private :: hdf_add_string
    procedure,private :: hdf_adda_string
  end type hdf5_file

contains
  !=============================================================================
  subroutine hdf_open_file(self,filename,status,action)
    !< Opens hdf5 file
    class(hdf5_file), intent(inout)    :: self
    character(*), intent(in)           :: filename
    character(*), intent(in), optional :: status
    character(*), intent(in), optional :: action

    character(:), allocatable :: lstatus, laction
    integer :: ierr
    logical :: is_hdf5

    !> Initialize FORTRAN interface.
    if(.not.self%is_init) call h5open_f(ierr)
    if (ierr /= 0 ) error stop 'Error: HDF5 library initialize Failed!'
    self%is_init = .TRUE.

    self%filename = filename
    inquire(file=filename, exist=is_hdf5)

    lstatus = 'old'
    if(present(status)) lstatus = to_lower(status)

    laction = 'rw'
    if(present(action)) laction = to_lower(action)

    select case(lstatus)
      case ('old')
        select case(laction)
          case('read','r')  !> Open an existing file.
            call h5fopen_f(filename,H5F_ACC_RDONLY_F,self%lid,ierr)
          case('write','readwrite','w','rw')
            if (is_hdf5) then
              call h5fopen_f(filename,H5F_ACC_RDWR_F,self%lid,ierr)
            else
              call h5fcreate_f(filename, H5F_ACC_TRUNC_F,self%lid,ierr)
            endif
          case default
            error stop 'Error: Unsupported action ->'// laction
          endselect
      case('new','replace')
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,self%lid,ierr)
      case default
        error stop 'Error: Unsupported status ->'// lstatus
    endselect

  end subroutine hdf_open_file
  !=============================================================================
  subroutine hdf_close_file(self, finalize)
    class(hdf5_file), intent(in) :: self
    logical, intent(in), optional:: finalize

    integer :: ierr

    !> close hdf5 file
    call h5fclose_f(self%lid, ierr)
    if(ierr /=0) error stop "Unable to close HDF5 file: "//self%filename

    !>  Close FORTRAN interface.
    if(present(finalize) .and. finalize) call h5close_f(ierr)
    if(ierr /=0) error stop "Unable to close HDF5 fortran interface!"
    ! this%is_init = .FALSE.

  end subroutine hdf_close_file
  !=============================================================================
  subroutine hdf_open_group(self,gname)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in)        :: gname

    integer :: ierr

    call h5gopen_f(self%lid, gname, self%gid, ierr)
    self%glid = self%lid
    self%lid  = self%gid

  end subroutine hdf_open_group
  !=============================================================================
  subroutine hdf_close_group(self)
    class(hdf5_file), intent(inout) :: self

    integer :: ierr

    call h5gclose_f(self%gid, ierr)
    self%lid = self%glid

  end subroutine hdf_close_group
  !=============================================================================
  function hdf_exist(self, name) result(exist)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: name
    logical :: exist

    integer :: ierr

    call h5lexists_f(self%lid, name, exist, ierr)

  end function hdf_exist
  !=============================================================================
  subroutine hdf_delete(self,name)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: name

    integer :: ierr

    call h5ldelete_f(self%lid, name, ierr)

  end subroutine hdf_delete
  !=============================================================================
  subroutine hdf_add_group(self, gname)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: gname    !< relative path to group

    integer(HID_T)  :: gid

    integer :: ierr, sp, ep, sl
    logical :: gexist

    sl = len(gname)
    sp = 1
    ep = 0

    do
       ep = index(gname(sp+1:sl), "/")

       ! no subgroup found
       if (ep == 0) exit

       ! check subgroup exists
       sp = sp + ep
       call h5lexists_f(self%lid, gname(1:sp-1), gexist, ierr)

       if(.not.gexist) then
         call h5gcreate_f(self%lid, gname(1:sp-1), gid, ierr)
         call h5gclose_f(gid, ierr)
       endif
    end do

  end subroutine hdf_add_group
  !=============================================================================
  subroutine hdf_adda_string(self, path, name, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: path
    character(*), intent(in) :: name
    character(*), intent(in) :: value

    integer :: ierr

    if(.not.self%exist(path)) call self%add(path)
    call h5ltset_attribute_string_f(self%lid, path, name, value, ierr)

  end subroutine hdf_adda_string
  !=============================================================================
  subroutine hdf_add_int(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in)      :: value

    integer(HID_T)  :: sid,did
    integer         :: ierr

    call self%add(dname)

    !> create dataspace
    call h5screate_f(H5S_SCALAR_F, sid, ierr)

    !> create dataset
    call h5dcreate_f(self%lid, dname, h5kind_to_type(kind(value),H5_INTEGER_KIND), sid, did, ierr)

    !> write dataset
    call h5dwrite_f(did, h5kind_to_type(kind(value),H5_INTEGER_KIND), value, int(shape(value),HSIZE_T), ierr)

    !> close space and dataset
    call h5dclose_f(did, ierr)

    call h5sclose_f(sid, ierr)


  end subroutine hdf_add_int
  !=============================================================================
  subroutine hdf_add_int1d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in)      :: value(:)

    integer         :: ierr

    call self%add(dname)


    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_INTEGER_KIND), value, ierr)

  end subroutine hdf_add_int1d
  !=============================================================================
  subroutine hdf_add_int2d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in)      :: value(:,:)

    integer         :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_INTEGER_KIND), value, ierr)

  end subroutine hdf_add_int2d
  !=============================================================================
  subroutine hdf_add_int3d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in)      :: value(:,:,:)

    integer         :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_INTEGER_KIND), value, ierr)

  end subroutine hdf_add_int3d
  !=============================================================================
  subroutine hdf_add_real(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(kind=dp), intent(in)      :: value

    integer(HID_T)  :: sid,did
    integer         :: ierr

    call self%add(dname)

    !> create dataspace
    call h5screate_f(H5S_SCALAR_F, sid, ierr)

    !> create dataset
    call h5dcreate_f(self%lid, dname, h5kind_to_type(kind(value),H5_REAL_KIND), sid, did, ierr)

    !> write dataset
    call h5dwrite_f(did, h5kind_to_type(kind(value),H5_REAL_KIND), value, int(shape(value),HSIZE_T), ierr)

    !> close space and dataset
    call h5dclose_f(did, ierr)

    call h5sclose_f(sid, ierr)


  end subroutine hdf_add_real
  !=============================================================================
  subroutine hdf_add_real1d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(kind=dp), intent(in) :: value(:)

    integer         :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_REAL_KIND), value, ierr)

  end subroutine hdf_add_real1d
  !=============================================================================
  subroutine hdf_add_real2d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(kind=dp), intent(in)      :: value(:,:)

    integer         :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_REAL_KIND), value, ierr)

  end subroutine hdf_add_real2d
  !=============================================================================
  subroutine hdf_add_real3d(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(kind=dp), intent(in)      :: value(:,:,:)

    integer         :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
      rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_REAL_KIND), value, ierr)

  end subroutine hdf_add_real3d
  !=============================================================================
  subroutine hdf_add_string(self,dname,value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in)     :: dname
    character(*), intent(in)     :: value

    integer         :: ierr

    call self%add(dname)
    call h5ltmake_dataset_string_f(self%lid, dname, value, ierr)

  end subroutine hdf_add_string
  !=============================================================================
  subroutine hdf_get_string(self,dname,value)
    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    character(:), intent(out), allocatable :: value

    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(character(dsize-1) :: value) !< string length = dsize -1
    call h5ltread_dataset_string_f(self%lid, dname, value, ierr)

  end subroutine hdf_get_string
  !=============================================================================
  subroutine hdf_get_int(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    integer, intent(out)             :: value

    integer(HID_T)  :: set_id
    integer :: ierr

    ! open dataset
    call h5dopen_f(self%lid, dname, set_id, ierr)

    ! read dataset
    call h5dread_f(set_id, h5kind_to_type(kind(value),H5_INTEGER_KIND), value,int(shape(value),HSIZE_T), ierr)

    ! close dataset
    call h5dclose_f(set_id, ierr)

  end subroutine hdf_get_int
  !=============================================================================
  subroutine hdf_get_int1d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    integer, intent(out),allocatable :: value(:)

    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_INTEGER_KIND), value, dims,  ierr)

  end subroutine hdf_get_int1d
  !=============================================================================
  subroutine hdf_get_int2d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    integer, intent(out),allocatable :: value(:,:)

    integer(HSIZE_T) :: dims(2)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1),dims(2)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_INTEGER_KIND), value, dims,  ierr)

  end subroutine hdf_get_int2d
  !=============================================================================
  subroutine hdf_get_int3d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    integer, intent(out),allocatable :: value(:,:,:)

    integer(HSIZE_T) :: dims(3)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1),dims(2),dims(3)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_INTEGER_KIND), value, dims,  ierr)

  end subroutine hdf_get_int3d
  !=============================================================================
  subroutine hdf_get_real(self, dname, value)

    class(hdf5_file), intent(in)  :: self
    character(*), intent(in)      :: dname
    real(kind=dp), intent(out)             :: value

    integer(HID_T)  :: set_id
    integer :: ierr

    ! open dataset
    call h5dopen_f(self%lid, dname, set_id, ierr)

    ! read dataset
    call h5dread_f(set_id, h5kind_to_type(kind(value),H5_REAL_KIND),&
                 & value,int(shape(value),HSIZE_T), ierr)

    ! close dataset
    call h5dclose_f(set_id, ierr)

  end subroutine hdf_get_real
  !=============================================================================
  subroutine hdf_get_real1d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    real(kind=dp), intent(out),allocatable :: value(:)

    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_REAL_KIND), value, dims,  ierr)

  end subroutine hdf_get_real1d
  !=============================================================================
  subroutine hdf_get_real2d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    real(kind=dp), intent(out),allocatable :: value(:,:)

    integer(HSIZE_T) :: dims(2)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1),dims(2)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_REAL_KIND), value, dims,  ierr)

  end subroutine hdf_get_real2d
  !=============================================================================
  subroutine hdf_get_real3d(self, dname, value)

    class(hdf5_file), intent(in)     :: self
    character(*), intent(in)         :: dname
    real(kind=dp), intent(out),allocatable :: value(:,:,:)

    integer(HSIZE_T) :: dims(3)
    integer(SIZE_T)  :: dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)

    allocate(value(dims(1),dims(2),dims(3)))

    call h5ltread_dataset_f(self%lid, dname, &
         & h5kind_to_type(kind(value),H5_REAL_KIND), value, dims,  ierr)

  end subroutine hdf_get_real3d

!----- Helper functions

  elemental function to_lower(str)
  ! can be trivially extended to non-ASCII
    character(*), intent(in) :: str
    character(len(str)) :: to_lower
    character(*), parameter :: lower="abcdefghijklmnopqrstuvwxyz", &
                               upper="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    integer :: i,j

    to_lower = str

    do concurrent (i = 1:len(str))
      j = index(upper,str(i:i))
      if (j > 0) to_lower(i:i) = lower(j:j)
    end do

  end function to_lower
  !=============================================================================

  subroutine h5read_f_1d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), allocatable, intent(out) :: value(:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%get(dname, value)
    call hdf%close()
    
  end subroutine h5read_f_1d

  subroutine h5read_f_2d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), allocatable, intent(out) :: value(:,:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%get(dname, value)
    call hdf%close()
    
  end subroutine h5read_f_2d

  subroutine h5read_f_3d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), allocatable, intent(out) :: value(:,:,:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%get(dname, value)
    call hdf%close()
    
  end subroutine h5read_f_3d

  subroutine h5write_f_1d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), intent(in) :: value(:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%add(dname, value)
    call hdf%close()
    
  end subroutine h5write_f_1d

  subroutine h5write_f_2d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), intent(in) :: value(:,:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%add(dname, value)
    call hdf%close()
    
  end subroutine h5write_f_2d

  subroutine h5write_f_3d(fname, dname, value)
    character(*), intent(in) :: fname
    character(*), intent(in) :: dname
    real(kind=dp), intent(in) :: value(:,:,:)

    type(hdf5_file) :: hdf
    integer :: ierr

    call hdf%open(fname)
    call hdf%add(dname, value)
    call hdf%close()
    
  end subroutine h5write_f_3d

end module hdf5_interface
