!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                 (there are currently many more authors!)
!                           (c) October 2023
!   
!     Changing History: Oct 2023, Initialize Codes
!
!=====================================================================

module measadj
  use my_mpi
  use para, ap => att_para_global
  use src_rec
  use model, am => att_model_global
  use para
  use model

  implicit none

  type, public :: att_measadj
    real(kind=dp)                                          :: srcx, srcy, chi
    integer,public                                         :: pidx, nrec
    real(kind=dp),dimension(:), allocatable                :: tt, vel, dist, recx, recy, tekio
    real(kind=dp),dimension(:), allocatable                :: weight
    real(kind=dp), dimension(:,:),allocatable              :: timetable
    integer, dimension(:), allocatable                     :: idx_path
    type(srcrec), pointer                                  :: sr
    contains
    ! procedure :: init => initialize_measadj
    procedure :: run_forward, run_adjoint, to_table, get_recs,distory, run_adjoint_density
  end type att_measadj

  contains

  subroutine get_recs(this, this_sr, pidx, src_name)
    class(att_measadj), intent(inout) :: this
    type(srcrec), target, intent(in) :: this_sr
    integer, intent(in) :: pidx
    character(len=*), intent(in) :: src_name
    
    this%sr => this_sr
    this%pidx = pidx
    this%timetable = zeros(am%n_xyz(1),am%n_xyz(2))
    
    call this%sr%stations%get_sta_pos(src_name, this%srcx, this%srcy)
    call this%sr%get_evt_gather(src_name, this%sr%periods(pidx),&
                           this%tt,this%vel,this%dist,this%weight,&
                           this%recx, this%recy, this%idx_path)
    this%nrec = size(this%tt)
  end subroutine get_recs

  subroutine distory(this)
    class(att_measadj), intent(out) :: this
    
  end subroutine distory

  subroutine run_forward(this, vel2d, m11, m22, m12, ref_t)
    class(att_measadj), intent(inout) :: this
    real(kind=dp), dimension(:,:),intent(in) :: vel2d, m11, m12, m22,ref_t
    real(kind=dp), dimension(:,:), allocatable :: vtmp

    vtmp = 1./vel2d
    call FSM_UW_PS_lonlat_2d(am%xgrids, am%ygrids,am%n_xyz(1),am%n_xyz(2),&
                             m11, m22, -m12,&
                             this%timetable,vtmp,this%srcx,this%srcy,ref_t)
    this%tekio = interp2(am%xgrids, am%ygrids, this%timetable, this%recx, this%recy)
    this%chi = 0.5*sum(this%weight*(this%tekio - this%tt)**2)
  end subroutine run_forward

  subroutine run_adjoint(this, m11, m22, m12, adjtable)
    class(att_measadj), intent(inout) :: this
    real(kind=dp), dimension(am%n_xyz(1),am%n_xyz(2)),intent(in) :: m11, m12, m22
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: adjtable
    integer :: i

    adjtable = zeros(am%n_xyz(1),am%n_xyz(2))
    do i = 1, this%nrec
      if (abs(this%tekio(i) - this%tt(i))/this%tt(i)>0.5) this%weight(i) = 0.
    enddo
    call FSM_O1_JSE_lonlat_2d(am%xgrids, am%ygrids,am%n_xyz(1),am%n_xyz(2),&
                              m11, m22, -m12,&
                              this%timetable, adjtable,this%recx, this%recy,&
                              this%weight*(this%tekio - this%tt), this%nrec)
    call mask(am%xgrids,am%ygrids,am%n_xyz(1),am%n_xyz(2),adjtable,this%srcx,this%srcy)
  end subroutine run_adjoint

  subroutine run_adjoint_density(this, m11, m22, m12, density)
    class(att_measadj), intent(inout) :: this
    real(kind=dp), dimension(am%n_xyz(1),am%n_xyz(2)),intent(in) :: m11, m12, m22
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: density
    integer :: i

    density = zeros(am%n_xyz(1),am%n_xyz(2))
    call FSM_O1_JSE_lonlat_2d(am%xgrids, am%ygrids,am%n_xyz(1),am%n_xyz(2),&
                              m11, m22, -m12,&
                              this%timetable, density,this%recx, this%recy,&
                              -1.0*this%weight, this%nrec)
    call mask(am%xgrids,am%ygrids,am%n_xyz(1),am%n_xyz(2),density,this%srcx,this%srcy)
  end subroutine run_adjoint_density

  subroutine to_table(this, tt_fwd)
    class(att_measadj), intent(inout) :: this
    real(kind=dp), dimension(:), intent(inout) :: tt_fwd
    integer :: i
    
    do i = 1, this%nrec
      tt_fwd(this%idx_path(i)) = this%tekio(i)
    enddo
  end subroutine


end module
