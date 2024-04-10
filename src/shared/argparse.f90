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
module argparse
  use shared_par
  use my_mpi, only: finalize_mpi
  implicit none

contains
  subroutine argparse_tomo(fname, isfwd)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    logical, intent(out) :: isfwd
    character(len=MAX_STRING_LEN),intent(out) :: fname
    integer :: i,nargs,m, nopt=2

    isfwd = .false.
    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i)) 
    enddo
    if (myrank==0 .and. (nargs==0 .or. any(args == '-h'))) then
      write(*,*)'Usage: surfatt_tomo -i para_file [-f] [-h]'
      write(*,*)''
      write(*,*)'Adjoint-state travel time tomography for surface wave'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i para_file  Path to parameter file in yaml format'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -f            Forward simulate travel time for surface wave instead of inversion, defaults to False'
      write(*,*)' -h            Print help message'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
      elseif(arg(1:2) == '-f') then
        m = m+1
        isfwd = .true.
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments'
    endif
  end subroutine argparse_tomo

  subroutine argparse_cb_fwd(fname, ncb, pert_vel, hmarg, anom_size)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    character(len=MAX_STRING_LEN),intent(out) :: fname
    integer, dimension(3), intent(out) :: ncb
    real(kind=dp), intent(out) :: pert_vel, hmarg, anom_size
    integer :: i,nargs,m,ier, nopt=5

    pert_vel = 0.08
    hmarg = 0.
    anom_size = 0.
    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i)) 
    enddo
    if (myrank==0 .and. (nargs==0 .or. any(args == '-h'))) then
      write(*,*)'Usage: surfatt_cb_fwd -i para_file -n nx/ny/nz [-h] '// &
                '[-m margin_degree] [-p pert] [-s anom_size_km]'
      write(*,*)''
      write(*,*)'Create checkerboard and forward simulate travel time for surface wave'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i para_file        Path to parameter file in yaml format'
      write(*,*)' -n nx/ny/nz         Number of anomalies along X, Y and Z direction'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -h                  Print help message'
      write(*,*)' -m margin_degree        Margin area in degree between anomaly as boundary, defaults to 0'
      write(*,*)' -p pert_vel         Magnitude of velocity perturbations, defaults to 0.08'
      write(*,*)' -s anom_size_km     size of anomalies at the top in km, default to uniform anomaly size in Z direction'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
      elseif(arg(1:2) == '-n') then
        m = m+1
        call parse_3string(args(i+1), ncb)
      elseif(arg(1:2) == '-p') then
        m = m+1
        read(args(i+1),*,iostat=ier) pert_vel
        if(ier/=0)stop 'Cannot parse argument -p'
      elseif(arg(1:2) == '-m') then
        m = m+1
        read(args(i+1),*,iostat=ier) hmarg
        if(ier/=0) stop 'Cannot parse argument -m'
      elseif(arg(1:2) == '-s') then
        m = m+1
        read(args(i+1),*,iostat=ier) anom_size
        if(ier/=0) stop 'Cannot parse argument -s'          
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments'
    endif
  end subroutine argparse_cb_fwd

  subroutine argparse_rotate_src_rec(fname, angle, center, outfname)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    character(len=MAX_STRING_LEN),intent(out) :: fname, outfname
    real(kind=dp), intent(out) :: angle, center(2)
    integer :: i,nargs,m,ier,nopt=4

    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i))
    enddo
    if (myrank==0 .and. (nargs==0 .or. any(args == '-h'))) then
      write(*,*)'Usage: surfatt_rotate_src_rec -i src_rec_file -a angle -c clat/clon [-h] [-o out_src_rec_file]'
      write(*,*)''
      write(*,*)'Rotate source and receiver locations by a given angle (anti-clockwise)'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i src_rec_file      Path to src_rec file in csv format'
      write(*,*)' -a angle             Angle in degree to rotate source and receiver locations'
      write(*,*)' -c clat/clon         Center of rotation in latitude and longitude'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -h                   Print help message'
      write(*,*)' -o out_file          Output file name, defaults to src_rec_file with "_rot" appended'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
        outfname = trim(adjustl(fname)) // '_rot'
      elseif(arg(1:2) == '-a') then
        m = m+1
        read(args(i+1),*,iostat=ier) angle
        if(ier/=0)stop 'Cannot parse angle in real format'
      elseif(arg(1:2) == '-c') then
        m = m+1
        call parse_2string_dp(args(i+1), center)
      elseif(arg(1:2) == '-o') then
        m = m+1
        outfname = args(i+1)
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments'
    endif
    
  end subroutine argparse_rotate_src_rec

  subroutine argparse_rotate_topo(fname, xrange, yrange, angle, center, outfname)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    character(len=MAX_STRING_LEN),intent(out) :: fname, outfname
    real(kind=dp), intent(out) :: angle, center(2), xrange(2), yrange(2)
    integer :: i,nargs,m,ier,nopt=6

    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i))
    enddo
    if (myrank==0 .and. (nargs==0 .or. any(args == '-h'))) then
      write(*,*)'Usage: surfatt_rotate_topo -i topo_file -a angle -c clat/clon -o out_topo_file' // &
                ' -x xmin/xmax -y ymin/ymax [-h]'
      write(*,*)''
      write(*,*)'Rotate topography by a given angle (anti-clockwise)'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i topo_file         Path to topography file in netcdf format'
      write(*,*)' -a angle             Angle in degree to rotate topography'
      write(*,*)' -c clat/clon         Center of rotation in latitude and longitude'
      write(*,*)' -o out_file          Output file name'
      write(*,*)' -x xmin/xmax         Range of new x coordinates'
      write(*,*)' -y ymin/ymax         Range of new y coordinates'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -h                   Print help message'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
      elseif(arg(1:2) == '-a') then
        m = m+1
        read(args(i+1),*,iostat=ier) angle
        if(ier/=0)stop 'Cannot parse angle in real format'
      elseif(arg(1:2) == '-c') then
        m = m+1
        call parse_2string_dp(args(i+1), center)
      elseif(arg(1:2) == '-o') then
        m = m+1
        outfname = args(i+1)
      elseif(arg(1:2) == '-x') then
        m = m+1
        call parse_2string_dp(args(i+1), xrange)
      elseif(arg(1:2) == '-y') then
        m = m+1
        call parse_2string_dp(args(i+1), yrange)
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments' 
    endif

  end subroutine argparse_rotate_topo

  subroutine argparse_tomo2d(fname, isfwd, ncb, pert_vel, hmarg)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    logical, intent(out) :: isfwd
    character(len=MAX_STRING_LEN),intent(out) :: fname
    integer, dimension(2), intent(out) :: ncb
    real(kind=dp), intent(out) :: pert_vel, hmarg
    integer :: i,nargs,m, nopt=6, ier

    isfwd = .false.
    ncb = [0, 0]
    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i)) 
    enddo
    if (myrank==0 .and. (nargs==0 .or. any(args == '-h'))) then
      write(*,*)'Usage: surfatt_tomo -i para_file [-f] [-h] '// &
                '[-n nx/ny] [-p pert_vel] [-m margin_km]'
      write(*,*)''
      write(*,*)'Adjoint-state travel time tomography for surface wave'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i para_file        Path to parameter file in yaml format'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -f                  Forward simulate travel time for surface wave instead of inversion, defaults to False'
      write(*,*)' -m margin_km        Margin area in km between anomaly as boundary, defaults to 0km'
      write(*,*)' -n nx/ny            Number of anomalies along X, Y direction'
      write(*,*)' -p pert_vel         Magnitude of velocity perturbations, defaults to 0.08'
      write(*,*)' -h                  Print help message'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
      elseif(arg(1:2) == '-f') then
        m = m+1
        isfwd = .true.
      elseif(arg(1:2) == '-n') then
        m = m+1
        call parse_2string_i(args(i+1), ncb)
      elseif(arg(1:2) == '-p') then
        m = m+1
        read(args(i+1),*,iostat=ier) pert_vel
        if(ier/=0)stop 'Cannot parse argument -p'
      elseif(arg(1:2) == '-m') then
        m = m+1
        read(args(i+1),*,iostat=ier) hmarg
        if(ier/=0) stop 'Cannot parse argument -m'
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments'
    endif
  end subroutine argparse_tomo2d

  subroutine parse_3string(input, nums)
    character(len=*), intent(in) :: input
    integer, dimension(3), intent(out) :: nums
    integer :: slash1, slash2, ierr

    ! Find the positions of the slashes
    slash1 = index(input, "/")
    slash2 = index(input(slash1 + 1 :), "/") + slash1

    ! Read numbers from the string
    read(input(1 : slash1 - 1), *, iostat=ierr) nums(1)
    if(ierr/=0) stop 'Cannot parse argument -n'
    read(input(slash1 + 1 : slash2 - 1), *, iostat=ierr) nums(2)
    if(ierr/=0) stop 'Cannot parse argument -n'
    read(input(slash2 + 1 :), *, iostat=ierr) nums(3)
    if(ierr/=0) stop 'Cannot parse argument -n'

  end subroutine parse_3string

  subroutine parse_2string_dp(input, nums)
    character(len=*), intent(in) :: input
    real(kind=dp), dimension(2), intent(out) :: nums
    integer :: slash, ierr

    ! Find the positions of the slashes
    slash = index(input, "/")

    ! Read numbers from the string
    read(input(1 : slash - 1), *, IOSTAT=ierr) nums(1)
    if(ierr/=0) then
      write(*,*)'Cannot parse argument '//trim(input)
      stop
    endif
    read(input(slash + 1 :), *, iostat=ierr) nums(2)
    if(ierr/=0) then
      write(*,*) 'Cannot parse argument '//trim(input)
      stop 
    endif
  end subroutine parse_2string_dp

  subroutine parse_2string_i(input, nums)
    character(len=*), intent(in) :: input
    integer, dimension(2), intent(out) :: nums
    integer :: slash, ierr

    ! Find the positions of the slashes
    slash = index(input, "/")

    ! Read numbers from the string
    read(input(1 : slash - 1), *, IOSTAT=ierr) nums(1)
    if(ierr/=0) then
      write(*,*)'Cannot parse argument '//trim(input)
      stop
    endif
    read(input(slash + 1 :), *, iostat=ierr) nums(2)
    if(ierr/=0) then
      write(*,*) 'Cannot parse argument '//trim(input)
      stop
    endif
    
  end subroutine parse_2string_i
end module argparse
