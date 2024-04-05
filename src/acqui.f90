!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: March 2024, Initialize Codes
!
!=====================================================================
module acqui
  use utils
  use para
  use grid, ag_ph => att_grid_global_ph, ag_gr => att_grid_global_gr
  use src_rec, sr_ph => src_rec_global_ph, sr_gr => src_rec_global_gr
  use decomposer, amd => att_mesh_decomposer_global
  use setup_att_log
  use measadj
  use surfker

  implicit none

  type, public :: att_acqui
    integer                                                :: itype,nsrc,istart,iend,iter
    type(att_grid), pointer                                :: ag
    type(SrcRec), pointer                                  :: sr
    logical                                                :: is_fwd
    character(len=MAX_STRING_LEN)                          :: model_fname, module='ACQUI',&
                                                              final_fname, gr_name, message
    real(kind=dp), dimension(:,:,:,:), pointer             :: adj_s, adj_density
    real(kind=dp), dimension(:,:,:,:), allocatable         :: adj_s_local, adj_density_local
    real(kind=dp), dimension(:,:,:,:), pointer             :: sen_vsRc, sen_vpRc, sen_rhoRc
    real(kind=dp), dimension(:,:,:), pointer               :: ker_beta, ker_density
    integer, dimension(:,:), pointer                       :: isrcs
    real(kind=dp)                                          :: chi0
    contains
    procedure :: init => att_acqui_init, scatter => att_scatter_src_gather, post_proc => post_processing_for_kernel
    procedure :: allocate_shm_arrays, prepare_inv, regularize_ker_density, post_proc_eikokernel,&
                 forward_simulate, depthkernel, combine_kernels
  end type att_acqui

  type(att_acqui), target, public                          :: att_acqui_global_ph, att_acqui_global_gr
  integer :: win_adj_s,win_adj_density,win_isrcs,&
            win_sen_vsRc, win_sen_vpRc, win_sen_rhoRc,&
            win_ker_beta, win_ker_density
contains
  subroutine att_acqui_init(this, itype, is_fwd)
    class(att_acqui), intent(inout) :: this
    integer, intent(in) :: itype
    logical, optional, intent(in) :: is_fwd
    logical :: is_fwd_local

    if (present(is_fwd)) then
      this%is_fwd = is_fwd
    else
      this%is_fwd = .false.
    endif

    this%itype = itype
    this%gr_name = trim(ap%data%gr_name(itype))
    if (itype == 1) then
      this%ag => ag_ph
      this%sr => sr_ph
    elseif(itype == 2) then
      this%ag => ag_gr
      this%sr => sr_gr
    endif
    call this%scatter()
    call this%allocate_shm_arrays()
  end subroutine att_acqui_init

  subroutine att_scatter_src_gather(this)
    class(att_acqui), intent(inout) :: this
    integer, dimension(:,:), allocatable :: isrcs
    integer, dimension(:), allocatable :: iperiods
    integer :: i, j, np

    iperiods = zeros(this%sr%nperiod)
    isrcs = zeros(this%sr%npath, 2)
    this%nsrc = 0
    if (myrank == 0) then
      do j = 1, this%sr%stations%nsta
        if (any(this%sr%evtname==this%sr%stations%staname(j))) then  
          call this%sr%get_periods_by_src(this%sr%stations%staname(j), iperiods, np)
          do i = 1, np
            this%nsrc = this%nsrc+1        
            isrcs(this%nsrc, 1) = iperiods(i)
            isrcs(this%nsrc, 2) = j
          enddo
        endif
      enddo
      write(this%message,'(a,i0," ",a,a,i0,a)') 'Scatter ',this%nsrc,&
            trim(this%gr_name),' events to ',mysize," processors"
      call write_log(this%message,1,this%module)
    endif
    call synchronize_all()
    call bcast_all_singlei(this%nsrc)
    call prepare_shm_array_i_2d(this%isrcs, this%nsrc, 2, win_isrcs)
    if (myrank == 0) this%isrcs(1:this%nsrc, :) = isrcs(1:this%nsrc, :)
    call synchronize_all()
    call scatter_all_i(this%nsrc,mysize,myrank,this%istart,this%iend)
  end subroutine att_scatter_src_gather

  subroutine prepare_inv(this)
    class(att_acqui), intent(inout) :: this

    this%adj_s_local = 0._dp
    this%adj_density_local = 0._dp
    if (myrank == 0) then
      this%sr%tt_fwd = 0._dp
      this%adj_s = 0._dp
      this%adj_density = 0._dp
      this%sen_vsRc = 0._dp
      this%sen_vpRc = 0._dp
      this%sen_rhoRc = 0._dp
      this%ker_beta = 0._dp
      this%ker_density = 0._dp
    endif
    call synchronize_all()
  end subroutine prepare_inv

  subroutine forward_simulate(this, chi_global, istotable, isadj, verbose)
    !> Compute the eikonal kernel for seismic tomography.
    !! Inputs:
    !!   this: an object of type att_acqui
    class(att_acqui), intent(inout) :: this
    type(att_measadj) :: ma
    real(kind=dp) :: chi_local
    logical, intent(in) :: istotable, isadj
    logical, optional, intent(in) :: verbose
    logical :: verbose_local
    real(kind=dp), intent(out) :: chi_global 
    real(kind=dp), dimension(:), allocatable :: local_tt
    real(kind=dp), dimension(:,:), allocatable :: adj, kden
    character(len=MAX_STRING_LEN) :: fname
    integer :: iz,i

    if (present(verbose)) then
      verbose_local = verbose
    else
      verbose_local = .true.
    endif

    chi_local = 0
    chi_global = 0
    local_tt = zeros(this%sr%npath)
    ! Forward simulation for surface wave velocity
    if (verbose_local) call write_log("Calculating surface wave velocity from Vs model...",1,this%module)
    call this%ag%fwdsurf(am%vs3d_opt)
    if (verbose_local) call write_log('This is measuring misfit and computing adjoint field for '//&
                            trim(this%gr_name)//'...',1,this%module)
    if ((this%iend-this%istart)>=0) then
      do i = this%istart, this%iend
        write(this%message, '(i0,a,F0.4,a,a)') myrank,' period: ',&
              this%sr%periods(this%isrcs(i,1)),' src_name: ',&
              this%sr%stations%staname(this%isrcs(i, 2))
        if (verbose_local) call write_log(this%message,0,this%module)
        ! get receivers for this source
        call ma%get_recs(this%sr,this%isrcs(i,1),this%sr%stations%staname(this%isrcs(i, 2)))
        ! forward simulation
        call ma%run_forward(this%ag%svel(this%isrcs(i,1),:,:), this%ag%m11(this%isrcs(i, 1),:,:),&
                            this%ag%m22(this%isrcs(i,1),:,:), this%ag%m12(this%isrcs(i,1),:,:),&
                            this%ag%ref_t(this%isrcs(i,1),:,:))
        ! sum chi
        chi_local = chi_local + ma%chi
        ! save synthetic tt to table
        if (istotable) call ma%to_table(local_tt)
        if (isadj) then
          ! measure adjoint
          call ma%run_adjoint(this%ag%m11(this%isrcs(i, 1),:,:),this%ag%m22(this%isrcs(i,1),:,:),&
                              this%ag%m12(this%isrcs(i,1),:,:),adj)
          ! kernel density
          call ma%run_adjoint_density(this%ag%m11(this%isrcs(i, 1),:,:), this%ag%m22(this%isrcs(i,1),:,:),&
                                      this%ag%m12(this%isrcs(i,1),:,:),kden)
          ! post proc of eikonal kernel
          call this%post_proc_eikokernel(this%isrcs(i,1), adj, ma%timetable)
        endif
        ! distribute measadj
        call ma%distory()
      enddo
    endif
    call synchronize_all()
    if (istotable) call sum_all_1Darray_dp(local_tt, this%sr%tt_fwd, this%sr%npath)
    ! reduce chi
    call sum_all_dp(chi_local, chi_global)
    call bcast_all_singledp(chi_global)
    call synchronize_all()
  end subroutine forward_simulate

  subroutine depthkernel(this)
    ! FILEPATH: /Users/xumijian/Codes/SurfATT/src/tomo.f90
    !> Compute the depth kernel for seismic tomography.
    !! Inputs:
    !!   this: an object of type att_acqui
    !! Outputs:
    !!   sen_vsRc: a 4D array of real numbers representing the sensitivity kernel for Vs
    !!   sen_vpRc: a 4D array of real numbers representing the sensitivity kernel for Vp
    !!   sen_rhoRc: a 4D array of real numbers representing the sensitivity kernel for density
    class(att_acqui), intent(inout) :: this
    real(kind=dp), dimension(:,:,:,:), allocatable :: loc_sen_vsRc, loc_sen_vpRc, loc_sen_rhoRc

    call write_log('Computing depth kernel...',1,this%module)
    ! compoute depth kernel
    call depthkernel_decomp(am%vs3d_opt,amd%loc_ix_start,amd%loc_ix_end,&
                            amd%loc_iy_start,amd%loc_iy_end,&                      
                            ap%data%iwave,ap%data%igr,this%sr%periods,am%zgrids,&
                            loc_sen_vsRc,loc_sen_vpRc,loc_sen_rhoRc)
    if (ap%topo%is_consider_topo) then
      call correct_depth(loc_sen_vsRc, loc_sen_vpRc, loc_sen_rhoRc,&
                        this%ag%topo_angle, amd%loc_ix_start,amd%loc_ix_end,&
                        amd%loc_iy_start,amd%loc_iy_end, am%zgrids)
    endif
    call synchronize_all()
    ! gather depth kernel
    call amd%collect_sen(loc_sen_vsRc, loc_sen_vpRc, loc_sen_rhoRc,&
                         this%sen_vsRc, this%sen_vpRc, this%sen_rhoRc)
    call synchronize_all()
  end subroutine depthkernel

  subroutine combine_kernels(this)
    ! FILEPATH: /Users/xumijian/Codes/SurfATT/src/tomo.f90
    !> Combine horizental and depth kernels.
    !! Inputs:
    !!   this: an object of type att_acqui
    class(att_acqui), intent(inout) :: this
    integer :: ip, i

    call write_log('Combining eikonal and surface wave kernels...',1,this%module)
    call sum_all_1Darray_dp(this%adj_s_local, this%adj_s, this%sr%nperiod*am%n_xyz(1)*am%n_xyz(2)*am%n_xyz(3))
    call sum_all_1Darray_dp(this%adj_density_local, this%adj_density, this%sr%nperiod*am%n_xyz(1)*am%n_xyz(2)*am%n_xyz(3))
        if (myrank == 0) then
      do ip = 1, this%sr%nperiod
        this%ker_beta = this%ker_beta -this%adj_s(ip,:,:,:) * this%sen_vsRc(ip,:,:,:)
        this%ker_beta = this%ker_beta -this%adj_s(ip,:,:,:) * this%sen_vpRc(ip,:,:,:) &
                        * dalpha_dbeta(am%vs3d_opt) 
        this%ker_beta = this%ker_beta -this%adj_s(ip,:,:,:) * this%sen_rhoRc(ip,:,:,:) &
                        * drho_dalpha(am%vp3d_opt) * dalpha_dbeta(am%vs3d_opt)
        this%ker_density = this%ker_density - this%adj_density(ip,:,:,:) * this%sen_vsRc(ip,:,:,:) - &
                           this%adj_density(ip,:,:,:) * this%sen_vpRc(ip,:,:,:) * dalpha_dbeta(am%vs3d_opt) - &
                           this%adj_density(ip,:,:,:) * this%sen_rhoRc(ip,:,:,:) * &
                           drho_dalpha(am%vp3d_opt) * dalpha_dbeta(am%vs3d_opt)
      enddo
      this%ker_beta = this%ker_beta / this%chi0
    endif
    call synchronize_all()
  end subroutine combine_kernels

  subroutine post_processing_for_kernel(this)
    !> This subroutine performs tomography inversion using gradient descent method.
    !! FILEPATH: /Users/xumijian/Codes/SurfATT/src/tomo.f90
    !! @param[inout] this an object of type att_acqui
    class(att_acqui), intent(inout) :: this
    integer :: nxinv,nyinv,nzinv,nset
    real(kind=dp) :: step, gkmax
    real(kind=dp), dimension(:), allocatable :: gk,gk_precond
    real(kind=dp), dimension(:,:,:), allocatable :: update, precond
    character(len=MAX_STRING_LEN) :: fname

    if (myrank == 0) then
      nset = ap%inversion%ncomponents
      nxinv = ap%inversion%n_inv_grid(1)
      nyinv = ap%inversion%n_inv_grid(2)
      nzinv = ap%inversion%n_inv_grid(3)
      ! initial inversion grid kernel
      call write_log('This is post processing for '//trim(this%gr_name)//' kernel...',1,this%module)
      gk = zeros(nset*nxinv*nyinv*nzinv)
      gk_precond = zeros(nset*nxinv*nyinv*nzinv)
      update = zeros(am%n_xyz(1),am%n_xyz(2),am%n_xyz(3))
      ! to inversion grids
      if (ap%inversion%kdensity_coe > 0) then
        call this%regularize_ker_density(precond)
      else
        precond = ones(am%n_xyz(1),am%n_xyz(2),am%n_xyz(3))
      endif
      ! this%ker_beta = this%ker_beta*precond
      call inv_grid_iso(am%xinv,am%yinv,am%zinv, nxinv, nyinv, nzinv,&
                        nset,this%ker_beta,am%n_xyz(1),am%n_xyz(2),am%n_xyz(3),&
                        gk,am%xgrids,am%ygrids,am%zgrids)
      call inv_grid_iso(am%xinv,am%yinv,am%zinv, nxinv, nyinv, nzinv,&
                        nset,precond,am%n_xyz(1),am%n_xyz(2),am%n_xyz(3),&
                        gk_precond,am%xgrids,am%ygrids,am%zgrids)
      call inv2fwd_iso(am%xinv,am%yinv,am%zinv,nxinv,nyinv,nzinv,&
                       ap%inversion%ncomponents,am%n_xyz(1),am%n_xyz(2),am%n_xyz(3), &
                       gk*gk_precond,am%xgrids,am%ygrids,am%zgrids,update)
      update = update / nset
      gradient_s = gradient_s + update*ap%data%weights(this%itype)
    endif
    call synchronize_all()

  end subroutine post_processing_for_kernel

  subroutine regularize_ker_density(this, precond)
    class(att_acqui), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), allocatable, intent(out) :: precond

    precond = zeros(am%n_xyz(1), am%n_xyz(2), am%n_xyz(3))
    precond = abs(this%ker_density)/maxval(abs(this%ker_density))
    where(precond<1e-2)
      precond = 0
    elsewhere
      precond = 1/precond**ap%inversion%kdensity_coe
    endwhere
    ! call save_npy('precond.npy', precond)
    ! precond = precond/maxval(abs(precond))
  end subroutine regularize_ker_density

  subroutine post_proc_eikokernel(this, pidx, adjtable, kden)
    class(att_acqui), intent(inout) :: this
    integer, intent(inout) :: pidx
    real(kind=dp), dimension(:,:),allocatable :: Tx, Ty
    ! real(kind=dp), dimension(:,:,:,:), allocatable :: adj_s_local, adj_xi_local, adj_eta_local
    real(kind=dp),  dimension(:,:), intent(in) :: adjtable, kden
    real(kind=dp), dimension(am%n_xyz(1), am%n_xyz(2)) :: vtmp, lat_corr
    integer :: i, j, k

    vtmp = 1/this%ag%svel(pidx, :,:)
        do i = 1, am%n_xyz(3)
      this%adj_s_local(pidx, :,:,i) = this%adj_s_local(pidx, :,:,i)+adjtable * vtmp**3
      this%adj_density_local(pidx, :,:,i) = this%adj_density_local(pidx, :,:,i)+kden
    enddo
  end subroutine post_proc_eikokernel

  subroutine allocate_shm_arrays(this)
    class(att_acqui), intent(inout) :: this

    if (.not. this%is_fwd) then
      allocate(this%adj_s_local(this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3)))
      allocate(this%adj_density_local(this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3)))
      call prepare_shm_array_dp_4d(this%adj_s, this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_adj_s)
      call prepare_shm_array_dp_4d(this%adj_density, this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_adj_density)
      call prepare_shm_array_dp_4d(this%sen_vsRc, this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_sen_vsRc)
      call prepare_shm_array_dp_4d(this%sen_vpRc, this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_sen_vpRc)
      call prepare_shm_array_dp_4d(this%sen_rhoRc, this%sr%nperiod, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_sen_rhoRc)
      call prepare_shm_array_dp_3d(this%ker_beta, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_ker_beta)
      call prepare_shm_array_dp_3d(this%ker_density, am%n_xyz(1), am%n_xyz(2), am%n_xyz(3), win_ker_density)
    endif
    call synchronize_all()
  end subroutine allocate_shm_arrays
end module