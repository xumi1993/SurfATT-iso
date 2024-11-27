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
!                       Jan 2024, Sphere Geometry
!
!=====================================================================

module para

  use my_mpi
  use shared_par
  use stdlib_logger, only: debug_level, information_level
  implicit none
  
  ! ---------------------- sections of parameter ------------------
  type, public :: para_data
    character(len=MAX_STRING_LEN)                          :: src_rec_file_ph, src_rec_file_gr
    integer                                                :: iwave, igr
    character(MAX_STRING_LEN)                              :: wave_name
    character(MAX_STRING_LEN),dimension(2)                 :: gr_name = ['ph', 'gr']
    real(kind=dp), dimension(2)                            :: weights
    logical,dimension(2)                                   :: vel_type
  end type para_data

  type para_output
    character(len=MAX_STRING_LEN), public                  :: output_path='OUTPUT_FILES'
    integer, public                                        :: log_level=1, verbose_level=1
  end type para_output

  type para_domain
    real(kind=dp), dimension(2), public                    :: ref_pos, depth
    real(kind=dp), public                                  :: agl=0.
    integer, public                                        :: num_grid_margin
    real(kind=dp), dimension(3), public                    :: interval
  end type para_domain

  type para_topo
    logical, public                                        :: is_consider_topo
    character(len=MAX_STRING_LEN), public                  :: topo_file
    real(kind=dp), public                                  :: wavelen_factor
  end type para_topo

  type, public :: para_inversion
    integer                                                :: ncomponents, n_inv_grid(3), &
                                                              init_model_type,niter,max_sub_niter=10, &
                                                              optim_method=0
    character(len=MAX_STRING_LEN)                          :: init_model_path
    character(len=MAX_NAME_LEN), dimension(2)              :: optim_name = ['CG   ', 'LBFGS']
    real(kind=dp)                                          :: vel_range(2), step_length, min_derr=0.001,&
                                                              kdensity_coe=0.5, maxshrink=0.8, sigma_2d=0.
  end type para_inversion

  type, public :: att_para
    type(para_data)                                        :: data
    type(para_output)                                      :: output
    type(para_domain)                                      :: domain
    type(para_topo)                                        :: topo
    type(para_inversion)                                   :: inversion
    logical                                                :: is_joint_grph=.false.
    real(kind=dp), dimension(3)                            :: d_xyz
    contains
    procedure :: read => read_par_file
    procedure, private :: check_src_rec_file, force_weight
  end type att_para
  ! --------------------- end sections of parameter ----------------

  type(att_para), public                                   :: att_para_global
  contains

  subroutine read_par_file(this, fname)
    ! -------------------------------
    ! read parameters from yaml file
    ! ------------------------------
    use yaml, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar

    class(att_para), intent(inout) :: this
    character(len=MAX_STRING_LEN), intent(in) :: fname
    character(len=error_length) :: error
    class(type_node), pointer :: root
    class(type_dictionary), pointer :: data_sec, output, domain,topo,inversion
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    type (type_error), pointer :: io_err
    integer :: stat
    character(len=MAX_STRING_LEN) :: errmsg

    if (myrank == 0) then
      root => parse(fname, error = error)
      if (error/='') call exit_mpi(myrank, error)

      select type (root)
      class is (type_dictionary)
        ! get data section
        data_sec => root%get_dictionary('data',required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%data%src_rec_file_ph = data_sec%get_string('src_rec_file_ph', default='', error=io_err)
        this%data%src_rec_file_gr = data_sec%get_string('src_rec_file_gr', default='', error=io_err)
        ! if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%data%iwave = data_sec%get_integer('iwave', error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        if (this%data%iwave == RAYLEIGH) then
          this%data%wave_name = 'rayleigh'
        else
          this%data%wave_name = 'love'
        endif
        list => data_sec%get_list('vel_type',required=.true.,error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call read_logi_list(list,this%data%vel_type)
        call this%check_src_rec_file()
        list => data_sec%get_list('weights',required=.true.,error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call read_real_list(list, this%data%weights)
        call this%force_weight()


        ! read output section
        output => root%get_dictionary('output', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%output%output_path = output%get_string('output_path', error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call EXECUTE_COMMAND_LINE('mkdir -p '//trim(this%output%output_path),&
                                  exitstat=stat, cmdmsg=errmsg)
        if (stat /= 0) then
          call exit_MPI(myrank, errmsg)
        endif
        log_fname = trim(this%output%output_path)//'/'//trim(log_basename)
        loglevel = output%get_integer('log_level', error=io_err)
        if (loglevel == 0) then
          this%output%log_level = debug_level ! in stdlib_logger
        else
          this%output%log_level = information_level
        endif
        this%output%verbose_level = output%get_integer('verbose_level', error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))

        ! read domain section
        domain => root%get_dictionary('domain',required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        list => domain%get_list('depth',required=.true.,error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call read_real_list(list, this%domain%depth)
        list => domain%get_list('interval',required=.true.,error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call read_real_list(list, this%domain%interval)
        this%domain%num_grid_margin = domain%get_integer('num_grid_margin', default=0, error=io_err)

        ! read topo section
        topo => root%get_dictionary('topo',required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%topo%is_consider_topo = topo%get_logical('is_consider_topo',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%topo%topo_file = topo%get_string('topo_file',default='',error=io_err)
        this%topo%wavelen_factor = topo%get_real('wavelen_factor',error=io_err)

        ! read inversion section
        inversion => root%get_dictionary('inversion', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%ncomponents = inversion%get_integer('ncomponents',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%min_derr = inversion%get_real('min_derr',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        list => inversion%get_list('n_inv_grid', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        call read_i_list(list, this%inversion%n_inv_grid)
        this%inversion%init_model_type = inversion%get_integer('init_model_type',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%init_model_path = inversion%get_string('init_model_path',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        list => inversion%get_list('vel_range', required=.true., error=io_err)
        call read_real_list(list, this%inversion%vel_range)
        this%inversion%niter = inversion%get_integer('niter',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%step_length = inversion%get_real('step_length',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%optim_method = inversion%get_integer('optim_method',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        if (this%inversion%optim_method == 1 .and. this%output%verbose_level < 1) then
          call exit_mpi(myrank, 'optim_method=1 requires verbose_level >= 1')
        end if
        this%inversion%maxshrink = inversion%get_real('maxshrink',error=io_err)
        if (associated(io_err)) call exit_mpi(myrank, trim(io_err%message))
        this%inversion%max_sub_niter = inversion%get_integer('max_sub_niter',default=10,error=io_err)
        this%inversion%kdensity_coe = inversion%get_real('kdensity_coe',default=0.0,error=io_err)
        this%inversion%sigma_2d = inversion%get_real('sigma_2d',default=0.0,error=io_err)

      end select
      call root%finalize()
      deallocate(root)
    endif ! if (myrank == 0) then
    call synchronize_all()
  
  !---------- broadcast all parameters ------------
    ! broadcast data
    call bcast_all(loglevel)
    call bcast_all(this%data%src_rec_file_ph)
    call bcast_all(this%data%src_rec_file_gr)
    call bcast_all(this%data%iwave)
    call bcast_all(this%data%igr)
    call bcast_all(this%data%gr_name,2,MAX_STRING_LEN)
    call bcast_all(this%data%vel_type,2)
    ! broadcast output
    call bcast_all(this%output%output_path)
    call bcast_all(this%output%log_level)
    call bcast_all(this%output%verbose_level)
    ! broadcast domain
    call bcast_all(this%domain%ref_pos, 2)
    call bcast_all(this%domain%depth, 2)
    call bcast_all(this%domain%interval, 3)
    call bcast_all(this%domain%agl)
    call bcast_all(this%domain%num_grid_margin)

    ! broadcast topo
    call bcast_all(this%topo%is_consider_topo)
    call bcast_all(this%topo%topo_file)
    call bcast_all(this%topo%wavelen_factor)
    ! broadcast inversion
    call bcast_all(this%inversion%ncomponents)
    call bcast_all(this%inversion%min_derr)
    call bcast_all(this%inversion%n_inv_grid, 3)
    call bcast_all(this%inversion%init_model_type)
    call bcast_all(this%inversion%init_model_path)
    call bcast_all(this%inversion%vel_range, 2)
    call bcast_all(this%inversion%niter)
    call bcast_all(this%inversion%step_length)
    call bcast_all(this%inversion%maxshrink)
    call bcast_all(this%inversion%kdensity_coe)
    call bcast_all(this%inversion%max_sub_niter)
    call bcast_all(this%inversion%optim_method)
    call bcast_all(this%inversion%sigma_2d)
    
    call synchronize_all()

  end subroutine read_par_file

  subroutine read_real_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    real(kind=dp), dimension(:), intent(out) :: list_out
    integer :: i
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_real(0.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_real_list

  subroutine read_logi_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    logical, dimension(:), intent(out) :: list_out
    integer :: i
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_logical(.true.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_logi_list

  subroutine read_i_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    integer, dimension(:), intent(out) :: list_out
    integer :: i
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_integer(0)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_i_list

  subroutine check_src_rec_file(this)
    class(att_para), intent(inout) :: this
    logical :: ioexist

    if(this%data%vel_type(1)) then
      inquire(file=this%data%src_rec_file_ph, exist=ioexist)
      if (.not. ioexist) call exit_MPI(myrank, 'No such file of '&
          //trim(this%data%src_rec_file_ph))
    endif
    if(this%data%vel_type(2)) then
      inquire(file=this%data%src_rec_file_gr, exist=ioexist)
      if (.not. ioexist) call exit_MPI(myrank, 'No such file of '&
          //trim(this%data%src_rec_file_gr))
    endif

  end subroutine check_src_rec_file

  subroutine force_weight(this)
    class(att_para), intent(inout) :: this
    
    if (abs(this%data%weights(1) + this%data%weights(2)-1.0) > 0.0001) then
      call exit_mpi(myrank, 'Summed weights of ph and gr should be 1.0')
    endif
    if(this%data%vel_type(1) .and. (.not. this%data%vel_type(2))) then
      this%data%weights(1) = 1.
      this%data%weights(2) = 0.
    elseif(this%data%vel_type(2) .and. (.not. this%data%vel_type(1))) then
      this%data%weights(1) = 0.
      this%data%weights(2) = 1.
    endif
  end subroutine force_weight
end module para
