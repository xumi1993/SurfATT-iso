module src_rec_parser
  use utils
  use csv_module
  use constants

  implicit none

  type SrcRecRaw
    real(kind=dp), dimension(:), allocatable, public       :: stla, stlo, stel,&
                                                              evla, evlo, evel,& 
                                                              periods_all, weight
    real(kind=dp), dimension(:), allocatable, public       :: tt, vel, dist
    integer, public                                        :: npath, nfield
    character(len=MAX_NAME_LEN), dimension(:), allocatable, public   :: evtname, staname, header
    contains
    procedure, public :: read_raw_src_rec_file, write_raw_src_rec_file, col_index
    procedure, private :: remove_el_from_header, order_header
  end type SrcRecRaw
contains
  subroutine read_raw_src_rec_file(this, fname)
    class(SrcRecRaw), intent(inout) :: this
    character(len=MAX_STRING_LEN), intent(in) :: fname
    type(csv_file) :: fp
    logical :: lerr
    integer :: i
  
    call fp%read(fname,header_row=1,status_ok=lerr)
    if (.not. lerr) then
      write(*,*) "No such file of "//trim(fname)
      stop 
    endif
    
    call fp%get_header(this%header,lerr)
    this%nfield = size(this%header, 1)
    call fp%get(this%col_index('tt'), this%tt, lerr)
    this%npath = size(this%tt)
    call fp%get(this%col_index('staname'), this%staname, lerr)
    call fp%get(this%col_index('stla'), this%stla, lerr)
    call fp%get(this%col_index('stlo'), this%stlo, lerr)
    call fp%get(this%col_index('evtname'), this%evtname, lerr)
    call fp%get(this%col_index('evla'), this%evla, lerr)
    call fp%get(this%col_index('evlo'), this%evlo, lerr)
    call fp%get(this%col_index('period'), this%periods_all, lerr)
    if (any('weight'==this%header)) then
      call fp%get(this%col_index('weight'), this%weight, lerr)
    else
      this%weight = ones(this%npath)
    endif
    if (any('dist'==this%header)) then
      call fp%get(this%col_index('dist'), this%dist, lerr)
    else
      allocate(this%dist(this%npath))
      do i = 1, this%npath
        this%dist(i) = gps2dist(this%stla(i),this%stlo(i),this%evla(i),this%evlo(i))
      enddo
    endif
    if (any('vel'==this%header)) then
      call fp%get(this%col_index('vel'), this%vel, lerr)
    else
      this%vel = this%dist/this%tt
    endif
    call this%remove_el_from_header()
    call this%order_header()
  end subroutine read_raw_src_rec_file

  subroutine write_raw_src_rec_file(this, fname)
    class(SrcRecRaw), intent(in) :: this
    character(len=*), intent(in) :: fname
    type(csv_file) :: fp
    logical :: status_ok
    integer :: i, n

    call fp%initialize(enclose_strings_in_quotes=.false.)
    call fp%open(fname,n_cols=this%nfield,status_ok=status_ok)
    if ( .not. status_ok) then
      write(*,*) 'Cannot open '//trim(fname)
      stop
    endif
    call fp%add(this%header, trim_str=.true.)
    call fp%next_row()
    do i = 1, this%npath
      call fp%add(this%tt(i), real_fmt='(F11.6)')
      call fp%add(this%staname(i), trim_str=.true.)
      call fp%add(this%stla(i), real_fmt='(F11.6)')
      call fp%add(this%stlo(i), real_fmt='(F11.6)')
      call fp%add(this%evtname(i), trim_str=.true.)
      call fp%add(this%evla(i), real_fmt='(F11.6)')
      call fp%add(this%evlo(i), real_fmt='(F11.6)')
      call fp%add(this%periods_all(i), real_fmt='(F0.4)')
      if (any('weight'==this%header)) call fp%add(this%weight(i), real_fmt='(F6.4)')
      if (any('dist'==this%header)) call fp%add(this%dist(i), real_fmt='(F11.6)')
      if (any('vel'==this%header)) call fp%add(this%vel(i), real_fmt='(F11.6)')
      call fp%next_row()
    enddo
    call fp%close(status_ok)
  end subroutine write_raw_src_rec_file

  subroutine order_header(this)
    class(SrcRecRaw), intent(inout) :: this

    integer :: i, j, n
    character(len=MAX_NAME_LEN), dimension(:), allocatable :: new_header

    n = 0
    allocate(new_header(this%nfield))
    do i = 1, this%nfield
      if (this%header(i) == 'tt') then
        new_header(1) = 'tt'
      elseif (this%header(i) == 'staname') then
        new_header(2) = 'staname'
      elseif (this%header(i) == 'stla') then
        new_header(3) = 'stla'
      elseif (this%header(i) == 'stlo') then  
        new_header(4) = 'stlo'
      elseif (this%header(i) == 'evtname') then
        new_header(5) = 'evtname'
      elseif (this%header(i) == 'evla') then
        new_header(6) = 'evla'
      elseif (this%header(i) == 'evlo') then
        new_header(7) = 'evlo'
      elseif (this%header(i) == 'period') then
        new_header(8) = 'period'
      elseif (this%header(i) == 'weight') then
        new_header(9) = 'weight'
      elseif (this%header(i) == 'dist') then
        new_header(10) = 'dist'
      elseif (this%header(i) == 'vel') then
        new_header(11) = 'vel'
      endif
    enddo
    this%header = new_header
  end subroutine order_header

  subroutine remove_el_from_header(this)
    class(SrcRecRaw), intent(inout) :: this
    integer :: i, j, n
    character(len=MAX_NAME_LEN), dimension(:), allocatable :: new_header

    n = 0
    do i = 1, this%nfield
      if (this%header(i) /= 'stel' .and. this%header(i) /= 'evel') then
        n = n + 1
      endif
    enddo
    allocate(new_header(n))
    j = 0
    do i = 1, this%nfield
      if (this%header(i) /= 'stel' .and. this%header(i) /= 'evel') then
        j = j + 1
        new_header(j) = this%header(i)
      endif
    enddo
    this%header = new_header
    this%nfield = n
    
  end subroutine remove_el_from_header

  function col_index(this, colname) result(idx)
    class(SrcRecRaw), intent(in) :: this
    character(len=*) :: colname
    integer :: idx, i
    
    ! find the index of the column in a list with the given name
    if (.not. any(colname == this%header)) then
      print *, "Column ", colname, " not found in the header"
      stop
    endif
    do i = 1, this%nfield
      if (trim(colname) == trim(this%header(i))) then
        idx = i
      endif
    enddo
  end function col_index
end module src_rec_parser