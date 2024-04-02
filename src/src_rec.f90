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

module src_rec
  use csv_module
  use constants
  use utils
  use para, ap => att_para_global
  use my_mpi
  use src_rec_parser
  use setup_att_log
  implicit none

  type, public :: Stations
    character(len=MAX_NAME_LEN),dimension(:), allocatable  :: staname
    real(kind=dp), dimension(:), allocatable               :: stla, stlo
    real(kind=dp), dimension(:), allocatable               :: stx, sty
    integer                                                :: nsta
    contains
    procedure :: get_sta_pos
  end type Stations

  type, public :: SrcRec
    real(kind=dp), dimension(:), pointer                  :: stla, stlo,&
                                                             evla, evlo,& 
                                                             periods_all, weight,&
                                                             tt, vel, dist
    integer                                                :: npath, nfield
    character(len=MAX_NAME_LEN), dimension(:), pointer     :: evtname, staname, header
    real(kind=dp), dimension(:), pointer                   :: periods, meanvel
    real(kind=dp), dimension(:), pointer                   :: tt_fwd
    integer                                                :: nperiod
    type(stations)                                         :: stations
    character(len=MAX_STRING_LEN)                          :: module = 'SRCREC', message
    contains
    procedure :: read => read_src_rec_file, to_csv             
    procedure :: get_sta, get_periods, get_mean_vel, get_evt_gather,get_periods_by_src                                           
  end type SrcRec

  type(srcrec), target, public                             :: src_rec_global_ph,src_rec_global_gr
  type(stations), target, public                           :: stations_global
  integer :: win_tt_fwd, win_vel, win_dist, win_weight, win_stla, win_stlo, win_stel,&
             win_evla, win_evlo, win_evel, win_periods_all, win_tt, &
             win_staname, win_evtname, win_periods_sr, win_meanvel, win_stx, win_sty

contains
  subroutine read_src_rec_file(this, type)
    class(SrcRec), intent(inout) :: this
    type(SrcRecRaw) :: srr
    character(len=MAX_STRING_LEN) :: fname
    integer, optional :: type

    if (present(type)) then
      if (type==0) then
        fname = ap%data%src_rec_file_ph
      else
        fname = ap%data%src_rec_file_gr
      endif
    else
      fname = ap%data%src_rec_file_ph
    endif

    if (myrank == 0) then
      call srr%read_raw_src_rec_file(fname)
      this%npath = srr%npath
      this%nfield = srr%nfield
    endif
    call synchronize_all()
    call bcast_all_singlei(this%npath)
    call bcast_all_singlei(this%nfield)
    
    ! allocate memory
    call prepare_shm_array_dp_1d(this%tt_fwd, this%npath, win_tt_fwd)
    call prepare_shm_array_dp_1d(this%vel, this%npath, win_vel)
    call prepare_shm_array_dp_1d(this%stla, this%npath, win_stla)
    call prepare_shm_array_dp_1d(this%stlo, this%npath, win_stlo)
    call prepare_shm_array_dp_1d(this%evla, this%npath, win_evla)
    call prepare_shm_array_dp_1d(this%evlo, this%npath, win_evlo)
    call prepare_shm_array_dp_1d(this%periods_all, this%npath, win_periods_all)
    call prepare_shm_array_dp_1d(this%weight, this%npath, win_weight)
    call prepare_shm_array_dp_1d(this%dist, this%npath, win_dist)
    call prepare_shm_array_dp_1d(this%tt, this%npath, win_tt)
    call prepare_shm_array_ch_1d(this%staname, this%npath, MAX_NAME_LEN, win_staname)
    call prepare_shm_array_ch_1d(this%evtname, this%npath, MAX_NAME_LEN, win_evtname)
    call prepare_shm_array_ch_1d(this%header, this%nfield, MAX_NAME_LEN, win_staname)

    if (myrank == 0) then
      this%vel = srr%vel
      this%stla = srr%stla
      this%stlo = srr%stlo
      this%evla = srr%evla
      this%evlo = srr%evlo
      this%periods_all = srr%periods_all
      this%weight = srr%weight
      this%dist = srr%dist
      this%tt = srr%tt
      this%staname = srr%staname
      this%evtname = srr%evtname
      this%header = srr%header
    endif
    call synchronize_all()

    ! get periods array
    call this%get_periods()

    ! get mean velocity
    call this%get_mean_vel()

    ! get station list
    call this%get_sta()

  end subroutine

  subroutine get_periods(this)
    use stdlib_sorting, only: sort
  
    class(SrcRec), intent(inout) :: this
    real(kind=dp), dimension(this%npath) :: temp
    integer :: i, count

    count = 1
    temp(1) = this%periods_all(1)
    if (myrank == 0) then
      do i = 2, this%npath
        if (.not. any(temp == this%periods_all(i))) then
          count = count + 1
          temp(count) = this%periods_all(i)
        endif
      enddo
      this%nperiod = count
      call sort(temp(1:count))
    endif
    call synchronize_all()
    call bcast_all_singlei(this%nperiod)
    call prepare_shm_array_dp_1d(this%periods, this%nperiod, win_periods_sr)
    if (myrank == 0) this%periods = temp(1:count)
    call synchronize_all()

  end subroutine get_periods

  subroutine get_periods_by_src(this, src_name, ipx, n)
    class(SrcRec), intent(inout) :: this
    character(len=MAX_NAME_LEN), intent(in) :: src_name
    real(kind=dp), dimension(this%nperiod) :: sta_periods
    integer, dimension(this%nperiod), intent(out) :: ipx
    integer, intent(out) :: n
    integer :: i

    sta_periods = 0
    n = 0
    do i = 1, this%npath
      if (this%evtname(i) == src_name) then
        if (any(sta_periods == this%periods_all(i))) cycle
        n = n + 1
        sta_periods(n) = this%periods_all(i)
        ipx(n) = findloc(this%periods,this%periods_all(i),1)
      endif
    enddo
  end subroutine get_periods_by_src

  subroutine get_sta(this)
    class(SrcRec), intent(inout) :: this
    character(len=MAX_NAME_LEN), dimension(this%npath*2) :: temp
    real(kind=dp), dimension(this%npath*2) :: templa, templo
    integer :: i, count
  
    ! get unique stations
    count = 1
    temp(1) = this%staname(1)
    templa(1) = this%stla(1)
    templo(1) = this%stlo(1)
    if (myrank == 0) then
      do i = 2, this%npath
        if (.not. any(temp == this%staname(i))) then
          count = count + 1
          temp(count) = this%staname(i)
          templa(count) = this%stla(i)
          templo(count) = this%stlo(i)
        endif
      enddo
      do i = 1, this%npath
        if (.not. any(temp == this%evtname(i))) then
          count = count + 1
          temp(count) = this%evtname(i)
          templa(count) = this%evla(i)
          templo(count) = this%evlo(i)
        endif
      enddo
      allocate(this%stations%staname(count),&
              this%stations%stla(count),&
              this%stations%stlo(count))
      this%stations%staname = temp(1:count)
      this%stations%stla = templa(1:count)
      this%stations%stlo = templo(1:count)
      this%stations%nsta = count
    endif
    call synchronize_all()
    call bcast_all_singlei(this%stations%nsta)
    if (myrank > 0) then
      allocate(this%stations%staname(this%stations%nsta),&
               this%stations%stla(this%stations%nsta),&
               this%stations%stlo(this%stations%nsta))
    endif
    call bcast_all_ch_array(this%stations%staname, this%stations%nsta, MAX_NAME_LEN)
    call bcast_all_dp(this%stations%stla, this%stations%nsta)
    call bcast_all_dp(this%stations%stlo, this%stations%nsta)
    allocate(this%stations%stx(this%stations%nsta))
    allocate(this%stations%sty(this%stations%nsta))
    this%stations%stx = this%stations%stlo
    this%stations%sty = this%stations%stla
    call synchronize_all()

  end subroutine get_sta

  subroutine get_mean_vel(this)
    class(SrcRec), intent(inout) :: this
    integer :: i, n, j
    real(kind=dp) :: sumval

    if (myrank == 0) then
      allocate(this%meanvel(this%nperiod))
      do i = 1, this%nperiod
        sumval = 0.
        n = 0
        do j = 1, this%npath
          if(this%periods_all(j) == this%periods(i)) then
            sumval = sumval + this%vel(j)
            n = n + 1
          endif
        enddo 
        this%meanvel(i) = sumval/n
      enddo
    endif
    call synchronize_all()
    if (myrank > 0) allocate(this%meanvel(this%nperiod))
    call bcast_all_dp(this%meanvel, this%nperiod)
    call synchronize_all()
  end subroutine get_mean_vel

  subroutine get_evt_gather(this, evtname, period, tt, vel, dist, weight, stx, sty, ipath)
    class(SrcRec), intent(inout) :: this
    character(len=*), intent(in) :: evtname
    real(kind=dp), intent(in) :: period
    real(kind=dp), dimension(:), allocatable, intent(out) :: tt, vel, dist, stx, sty
    real(kind=dp), dimension(:), allocatable, intent(out) :: weight
    integer, dimension(:), allocatable, intent(out) :: ipath
    integer :: i, n, idx, indices(this%npath), j

    n = 0
    do i = 1, this%npath
      if (this%evtname(i) == evtname .and. this%periods_all(i)==period) then
        n = n + 1
        indices(n) = i
      endif
    enddo
    allocate(tt(n), vel(n), dist(n), stx(n), sty(n), weight(n),ipath(n))
    do i = 1, n
      idx = indices(i)
      ipath(i) = idx
      tt(i) = this%tt(idx)
      vel(i) = this%vel(idx)
      dist(i) = this%dist(idx)
      weight(i) = this%weight(idx)
      do j = 1, this%stations%nsta
        if (this%stations%staname(j) == this%staname(idx)) then
          stx(i) = this%stations%stx(j)
          sty(i) = this%stations%sty(j)
        endif
      enddo
    enddo
  end subroutine get_evt_gather

  subroutine to_csv(this, fname)
    class(srcrec), intent(in) :: this
    character(len=*), intent(in) :: fname
    type(csv_file) :: fp
    logical :: status_ok
    integer :: i, nfield=13,n

    if (myrank == 0) then
      call write_log('Writing data to: '//trim(fname),1,this%module)
      n = size(this%header)
      ! if (n==10) then
        ! nfield = n + 3
      ! else
        ! nfield = 13
      ! endif
      call fp%initialize(enclose_strings_in_quotes=.false.)
      call fp%open(fname,n_cols=nfield,status_ok=status_ok)
      if ( .not. status_ok) call exit_MPI(myrank, 'Cannot open '//trim(fname))
      call fp%add(this%header, trim_str=.true.)
      if (.not. any(this%header == 'weight')) call fp%add('weight', trim_str=.true.)
      if (.not. any(this%header == 'dist')) call fp%add('dist', trim_str=.true.)
      if (.not. any(this%header == 'vel')) call fp%add('vel', trim_str=.true.)
      call fp%next_row()
      do i = 1, this%npath
        call fp%add(this%tt_fwd(i), real_fmt='(F0.6)')
        call fp%add(this%staname(i), trim_str=.true.)
        call fp%add(this%stla(i), real_fmt='(F9.4)')
        call fp%add(this%stlo(i), real_fmt='(F9.4)')
        call fp%add(this%evtname(i), trim_str=.true.)
        call fp%add(this%evla(i), real_fmt='(F9.4)')
        call fp%add(this%evlo(i), real_fmt='(F9.4)')
        call fp%add(this%periods_all(i), real_fmt='(F0.4)')
        call fp%add(this%weight(i), real_fmt='(F6.4)')
        call fp%add(this%dist(i), real_fmt='(F12.4)')
        call fp%add(this%dist(i)/this%tt_fwd(i), real_fmt='(F10.6)')
        call fp%next_row()
      enddo
      call fp%close(status_ok)
    endif! if(myrank==0)
    call synchronize_all()
  end subroutine to_csv

  subroutine get_sta_pos(this, src_name, stax, stay)
    class(stations), intent(inout) :: this
    character(len=MAX_NAME_LEN), intent(in) :: src_name
    real(kind=dp), intent(out) :: stax, stay
    integer :: i

    do i = 1, this%nsta
      if (this%staname(i) == src_name) then
        stax = this%stx(i)
        stay = this%sty(i)
      endif
    enddo
  end subroutine get_sta_pos

  subroutine merge_sta()
    integer :: i, j, n, nph, ngr
    ! character(len=MAX_NAME_LEN), dimension(:), allocatable :: compsta
    ! real(kind=dp), dimension(:), allocatable :: stax, stay, stla, stlo, stel
    integer,dimension(:), allocatable :: idx

    if (myrank == 0) then
      nph = src_rec_global_ph%stations%nsta
      ngr = src_rec_global_gr%stations%nsta
      if (all(ap%data%vel_type)) then
        n = 0
        allocate(idx(nph))
        do i = 1, nph
          if (.not. any(src_rec_global_gr%stations%staname&
              ==src_rec_global_ph%stations%staname(i))) then
              n = n + 1
              idx(n) = i
          endif
        enddo
        if (n/=0) then
          stations_global%nsta = ngr+n
          allocate(stations_global%staname(ngr+n),stations_global%stx(ngr+n),&
                   stations_global%sty(ngr+n),stations_global%stla(ngr+n),&
                   stations_global%stlo(ngr+n))
          stations_global%staname(1:ngr) = src_rec_global_gr%stations%staname(1:ngr)
          stations_global%stx(1:ngr) = src_rec_global_gr%stations%stx(1:ngr)
          stations_global%sty(1:ngr) = src_rec_global_gr%stations%sty(1:ngr)
          stations_global%stla(1:ngr) = src_rec_global_gr%stations%stla(1:ngr)
          stations_global%stlo(1:ngr) = src_rec_global_gr%stations%stlo(1:ngr)
          ! stations_global%stel(1:ngr) = src_rec_global_gr%stations%stel(1:ngr)
          do i = 1, n
            stations_global%staname(ngr+i) = src_rec_global_ph%stations%staname(idx(i))
            stations_global%stx(ngr+i) = src_rec_global_ph%stations%stx(idx(i))
            stations_global%sty(ngr+i) = src_rec_global_ph%stations%sty(idx(i))
            stations_global%stla(ngr+i) = src_rec_global_ph%stations%stla(idx(i))
            stations_global%stlo(ngr+i) = src_rec_global_ph%stations%stlo(idx(i))
            ! stations_global%stel(ngr+i) = src_rec_global_ph%stations%stel(idx(i))
          enddo
        else
          stations_global%nsta = ngr
          allocate(stations_global%staname(ngr),stations_global%stx(ngr),&
                   stations_global%sty(ngr),stations_global%stla(ngr),&
                   stations_global%stlo(ngr))
          stations_global%staname(1:ngr) = src_rec_global_gr%stations%staname(1:ngr)
          stations_global%stx(1:ngr) = src_rec_global_gr%stations%stx(1:ngr)
          stations_global%sty(1:ngr) = src_rec_global_gr%stations%sty(1:ngr)
          stations_global%stla(1:ngr) = src_rec_global_gr%stations%stla(1:ngr)
          stations_global%stlo(1:ngr) = src_rec_global_gr%stations%stlo(1:ngr)
          ! stations_global%stel(1:ngr) = src_rec_global_gr%stations%stel(1:ngr)
        endif
      endif ! (all(ap%data%vel_type)) 
    endif ! myrank==0
    call synchronize_all()
    if (all(ap%data%vel_type)) then
      call bcast_all_singlei(stations_global%nsta)
      if (myrank > 0) then
        allocate(stations_global%staname(stations_global%nsta),&
                 stations_global%stx(stations_global%nsta),&
                 stations_global%sty(stations_global%nsta),&
                 stations_global%stla(stations_global%nsta),&
                 stations_global%stlo(stations_global%nsta))
      endif
      call bcast_all_ch_array(stations_global%staname,&
          stations_global%nsta,MAX_NAME_LEN)
      call bcast_all_dp(stations_global%stla, stations_global%nsta)
      call bcast_all_dp(stations_global%stlo, stations_global%nsta)
      call bcast_all_dp(stations_global%stx, stations_global%nsta)
      call bcast_all_dp(stations_global%sty, stations_global%nsta)
    elseif (ap%data%vel_type(1)) then
      stations_global = src_rec_global_ph%stations
    elseif (ap%data%vel_type(2)) then
      stations_global = src_rec_global_gr%stations
    endif
  end subroutine merge_sta
end module
