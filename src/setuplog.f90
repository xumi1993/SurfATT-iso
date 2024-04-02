module setup_att_log
  use para, ap => att_para_global
  use stdlib_logger, logger => global_logger
  implicit none

contains
  subroutine initialize_files()
    character(len=MAX_STRING_LEN) :: fname
    integer :: stat

    fname = trim(ap%output%output_path)//"/"//trim(modfile)
    open(unit=IOUT, iostat=stat, file=fname, status='old')
    if (stat == 0) close(IOUT, status='delete')
    fname = trim(ap%output%output_path)//'/final_model.h5'
    open(unit=IOUT, iostat=stat, file=fname, status='old')
    if (stat == 0) close(IOUT, status='delete')
    fname = trim(ap%output%output_path)//'/initial_model.h5'
    open(unit=IOUT, iostat=stat, file=fname, status='old')
    if (stat == 0) close(IOUT, status='delete')
    fname = trim(ap%output%output_path)//'/target_model.h5'
    open(unit=IOUT, iostat=stat, file=fname, status='old')
    if (stat == 0) close(IOUT, status='delete')
  end subroutine initialize_files

  subroutine setuplog()
    integer :: stat

    ! if (myrank == 0) then
    call logger%add_log_file(log_fname, LID, stat=stat)      
    call logger%configure(level=ap%output%log_level, time_stamp=.true.)
    ! end
    call bcast_all_singlei(LID)
    
  end subroutine setuplog

  subroutine write_log(msg, level, module_name)
    character(len=*), intent(in) :: msg, module_name
    integer :: level
    character(len=5) :: level_msg
    character(2) :: spliter=': '
    character(23) :: time_stamp
    character(8)  :: date
    character(10) :: time

    call date_and_time( date, time )

    time_stamp(1:4)   = date(1:4)
    time_stamp(5:5)   = '-'
    time_stamp(6:7)   = date(5:6)
    time_stamp(8:8)   = '-'
    time_stamp(9:10)  = date(7:8)
    time_stamp(11:11) = ' '
    time_stamp(12:13) = time(1:2)
    time_stamp(14:14) = ':'
    time_stamp(15:16) = time(3:4)
    time_stamp(17:17) = ':'
    time_stamp(18:23) = time(5:10)

    if (level == 0) then
      level_msg = 'DEBUG'
    elseif(level == 1) then
      level_msg = 'INFO'
    elseif(level == 2) then
      level_msg = 'WARN'
    elseif(level == 3) then
      level_msg = 'ERROR'
    endif

    if (myrank == 0 .and. level>0) then
      call logger % log_message( msg,                  &
                                module = module_name,  &
                                prefix = trim(level_msg) )
      call flush(LID)
    endif
    
    if (level >= loglevel) then
      if (level > 0) then
        if (myrank==0) write(*, *) time_stamp//spliter//trim(module_name)//spliter//&
              trim(level_msg)//spliter//trim(msg)
      else
        write(*, *) time_stamp//spliter//trim(module_name)//spliter//&
              trim(level_msg)//spliter//trim(msg)
      endif
    endif

  end subroutine write_log
end module setup_att_log