module decomposer

  use shared_par
  use para, ap => att_para_global
  use model, am => att_model_global

  implicit none

  type, public :: att_decomposer
    integer                                                :: glob_px, glob_py, loc_nx, loc_ny
    integer                                                :: loc_ix_start, loc_ix_end, loc_iy_start, loc_iy_end
    integer, dimension(:,:), allocatable                   :: glob_ix, glob_iy
    contains
    procedure :: init => init_decomposer, decompose => partial_decompose_mesh, collect => collect_mesh
  end type att_decomposer

contains

  subroutine init_decomposer(this)
    class(att_decomposer), intent(inout) :: this
    integer :: f1, f2, istart, iend
    integer, dimension(mysize,2) :: loc_ix, loc_iy
    
    call close_factors(mysize, f1, f2)
    if (am%n_xyz(1) <= am%n_xyz(2)) then
      this%glob_px = f1
      this%glob_py = f2
    else
      this%glob_px = f2
      this%glob_py = f1
    end if

    allocate(this%glob_ix(mysize, 2))
    allocate(this%glob_iy(mysize, 2))
    this%loc_ix_start = 1 + (mod(myrank, this%glob_px)) * (am%n_xyz(1) / this%glob_px)
    this%loc_ix_end = this%loc_ix_start + (am%n_xyz(1) / this%glob_px) - 1
    this%loc_iy_start = 1 + (myrank) / this%glob_px * (am%n_xyz(2) / this%glob_py)
    this%loc_iy_end = this%loc_iy_start + (am%n_xyz(2) / this%glob_py) - 1
    loc_ix(myrank+1, 1) = this%loc_ix_start
    loc_ix(myrank+1, 2) = this%loc_ix_end
    loc_iy(myrank+1, 1) = this%loc_iy_start
    loc_iy(myrank+1, 2) = this%loc_iy_end
    call sum_all_1Darray_i(loc_ix, this%glob_ix, mysize*2)
    call sum_all_1Darray_i(loc_iy, this%glob_iy, mysize*2)

  end subroutine init_decomposer

  function partial_decompose_mesh(this, data) result(partial_data)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: data
    real(kind=dp), dimension(:,:,:), allocatable :: partial_data

    allocate(partial_data(this%loc_nx, this%loc_ny, am%n_xyz(3)))
    partial_data = data(this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :)

  end function partial_decompose_mesh

  function collect_mesh(this, data) result(gathered_data)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: data
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_data
    real(kind=dp), dimension(:,:,:), allocatable :: gathered_data
    integer :: i, ierr, nx, ny, tag=100

    if (myrank == 0) then
      allocate(gathered_data(am%n_xyz(1), am%n_xyz(2), am%n_xyz(3)))
      gathered_data(this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :) = data
      do i = 1, mysize - 1
        nx = this%glob_ix(i+1, 2) - this%glob_ix(i+1, 1) + 1
        ny = this%glob_iy(i+1, 2) - this%glob_iy(i+1, 1) + 1
        allocate(tmp_data(nx, ny, am%n_xyz(3)))
        call recv_dp(tmp_data, this%loc_nx * this%loc_ny * am%n_xyz(3), i, tag)
        gathered_data(this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2), :) = tmp_data
        deallocate(tmp_data)
      end do
    else
      allocate(tmp_data(this%loc_nx, this%loc_ny, am%n_xyz(3)))
      tmp_data = data
      call send_dp(tmp_data, this%loc_nx * this%loc_ny * am%n_xyz(3), 0, tag)
      deallocate(tmp_data)
    end if
    
  end function collect_mesh

  subroutine close_factors(num, f1, f2)
    integer :: num, f1, f2, i

    f1 = 1; f2 = num
    do i = 2, int(sqrt(real(num)))
      if (mod(num, i) == 0) then
        f1 = i
        f2 = num / i
        return
      end if
    end do
  end subroutine close_factors
end module