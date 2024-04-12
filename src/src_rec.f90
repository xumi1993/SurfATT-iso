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
    character(len=MAX_NAME_LEN),dimension(:), pointer  :: staname
    real(kind=dp), dimension(:), pointer               :: stla, stlo
    integer                                            :: nsta
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
  type(Stations), target, public                           :: stations_global
  integer :: win_tt_fwd, win_vel, win_dist, win_weight, win_stla, win_stlo,&
             win_evla, win_evlo, win_periods_all, win_tt, &
             win_staname, win_evtname, win_periods_sr, win_meanvel,&
             win_sta_name, win_sta_la, win_sta_lo, win_glob_sta,&
             win_glob_stla, win_glob_stlo

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

    call write_log('Reading source receiver file: '//trim(fname),1,this%module)

    if (local_rank == 0) then
      call srr%read_raw_src_rec_file(fname)
      this%npath = srr%npath
      this%nfield = srr%nfield
    endif
    call synchronize_all()
    call bcast_all(this%npath)
    call bcast_all(this%nfield)
    
    ! allocate memory
    ! this%tt_fwd = zeros(this%npath)
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

    if (local_rank == 0) then
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
    ! use stdlib_sorting, only: sort
  
    class(SrcRec), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: temp
    integer :: i, count

    count = 1
    if (local_rank == 0) then
      call append(temp, this%periods_all(1))
      do i = 2, this%npath
        if (.not. any(temp == this%periods_all(i))) then
          count = count + 1
          call append(temp, this%periods_all(i))
        endif
      enddo
      this%nperiod = count
      temp(1:count) = sort(temp(1:count))
    endif
    call synchronize_all()
    call bcast_all(this%nperiod)
    call prepare_shm_array_dp_1d(this%periods, this%nperiod, win_periods_sr)
    if (local_rank == 0) then
      this%periods = temp(1:count)
      if (this%periods(this%nperiod)*1.4 > ap%domain%depth(2)) then
        call write_log('the depth of the model smaller than 1.4 * maximum period',2,this%module)
      endif
      deallocate(temp)
    endif
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
    character(len=MAX_NAME_LEN), dimension(:), allocatable :: temp
    real(kind=dp), dimension(:), allocatable :: templa, templo
    integer :: i, count
  
    ! get unique stations
    if (myrank == 0) then
      count = 1
      call append(temp, this%staname(1))
      call append(templa, this%stla(1))
      call append(templo, this%stlo(1))
      do i = 2, this%npath
        if (.not. any(temp == this%staname(i))) then
          count = count + 1
          call append(temp, this%staname(i))
          call append(templa, this%stla(i))
          call append(templo, this%stlo(i))
        endif
      enddo
      do i = 1, this%npath
        if (.not. any(temp == this%evtname(i))) then
          count = count + 1
          call append(temp, this%evtname(i))
          call append(templa, this%evla(i))
          call append(templo, this%evlo(i))
        endif
      enddo
      this%stations%nsta = count
    endif
    call synchronize_all()
    call bcast_all(this%stations%nsta)
    call prepare_shm_array_ch_1d(this%stations%staname, this%stations%nsta, MAX_NAME_LEN, win_staname)
    call prepare_shm_array_dp_1d(this%stations%stla, this%stations%nsta, win_stla)
    call prepare_shm_array_dp_1d(this%stations%stlo, this%stations%nsta, win_stlo)
    if (local_rank == 0) then
      this%stations%staname = temp(1:count)
      this%stations%stla = templa(1:count)
      this%stations%stlo = templo(1:count)
    endif
    call synchronize_all()

  end subroutine get_sta

  subroutine get_mean_vel(this)
    class(SrcRec), intent(inout) :: this
    integer :: i, n, j
    real(kind=dp) :: sumval
    
    call prepare_shm_array_dp_1d(this%meanvel, this%nperiod, win_meanvel)
    if (local_rank == 0) then
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
  end subroutine get_mean_vel

  subroutine get_evt_gather(this, evtname, period, tt, vel, dist, weight, stlo, stla, ipath)
    class(SrcRec), intent(inout) :: this
    character(len=*), intent(in) :: evtname
    real(kind=dp), intent(in) :: period
    real(kind=dp), dimension(:), allocatable, intent(out) :: tt, vel, dist, stlo, stla
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
    allocate(tt(n), vel(n), dist(n), stlo(n), stla(n), weight(n),ipath(n))
    do i = 1, n
      idx = indices(i)
      ipath(i) = idx
      tt(i) = this%tt(idx)
      vel(i) = this%vel(idx)
      dist(i) = this%dist(idx)
      weight(i) = this%weight(idx)
      do j = 1, this%stations%nsta
        if (this%stations%staname(j) == this%staname(idx)) then
          stlo(i) = this%stations%stlo(j)
          stla(i) = this%stations%stla(j)
        endif
      enddo
    enddo
  end subroutine get_evt_gather

  subroutine to_csv(this, fname)
    class(srcrec), intent(in) :: this
    character(len=*), intent(in) :: fname
    type(csv_file) :: fp
    logical :: status_ok
    integer :: i, nfield=11,n

    if (myrank == 0) then
      call write_log('Writing data to: '//trim(fname),1,this%module)
      n = size(this%header)
      if (n==8) then
        nfield = n + 3
      else
        nfield = 11
      endif
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

    i = find_loc(this%staname, src_name)
    stay = this%stla(i)
    stax = this%stlo(i)
  end subroutine get_sta_pos

  subroutine merge_sta()
    integer :: i, j, nsta, nph, ngr
    character(len=MAX_NAME_LEN), dimension(:), allocatable :: tmpsta
    integer,dimension(:), allocatable :: idx
    real(kind=dp), dimension(:), allocatable :: tmpla, tmplo

    if (all(ap%data%vel_type)) then
      if (myrank == 0) then
        nph = src_rec_global_ph%stations%nsta
        ngr = src_rec_global_gr%stations%nsta
        tmpla = src_rec_global_ph%stations%stla
        tmplo = src_rec_global_ph%stations%stlo
        tmpsta = src_rec_global_ph%stations%staname
        do i = 1, ngr
          if (.not. any(tmpsta == src_rec_global_gr%stations%staname(i))) then
            call append(tmpsta, src_rec_global_gr%stations%staname(i))
            call append(tmpla, src_rec_global_gr%stations%stla(i))
            call append(tmplo, src_rec_global_gr%stations%stlo(i))
          endif
        enddo
        stations_global%nsta = size(tmpla)
      endif ! myrank==0
      call synchronize_all()
      call bcast_all(stations_global%nsta)
      call prepare_shm_array_dp_1d(stations_global%stla, stations_global%nsta, win_glob_stla)
      call prepare_shm_array_dp_1d(stations_global%stlo, stations_global%nsta, win_glob_stlo)
      call prepare_shm_array_ch_1d(stations_global%staname,stations_global%nsta,&
                                   MAX_NAME_LEN, win_glob_sta)
      if (myrank == 0) then
        stations_global%staname = tmpsta
        stations_global%stla = tmpla
        stations_global%stlo = tmplo
      endif
      call sync_from_main_rank(stations_global%staname, stations_global%nsta, MAX_NAME_LEN)
      call sync_from_main_rank(stations_global%stla, stations_global%nsta)
      call sync_from_main_rank(stations_global%stlo, stations_global%nsta)
    elseif (ap%data%vel_type(1)) then
      stations_global = src_rec_global_ph%stations
    elseif (ap%data%vel_type(2)) then
      stations_global = src_rec_global_gr%stations
    endif !     
    call synchronize_all()
  end subroutine merge_sta
end module
