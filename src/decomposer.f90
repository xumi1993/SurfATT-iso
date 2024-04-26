!=====================================================================
!
!                           S u r f A T T
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu @
!                     Nanyang Technological University
!                           (c) October 2023
!   
!     Changing History: Apl 2024, Initialize Codes
!
!=====================================================================
module decomposer

  use shared_par
  use utils
  use para, ap => att_para_global

  implicit none

  type, public :: att_decomposer
    integer                                                :: glob_px, glob_py, loc_nx, loc_ny
    integer                                                :: loc_ix_start, loc_ix_end, loc_iy_start, loc_iy_end
    integer, dimension(:,:), allocatable                   :: glob_ix, glob_iy
    contains
    procedure :: init => init_decomposer, decompose => partial_decompose_mesh
    procedure :: collect_mesh, collect_grid, collect_sen
  end type att_decomposer

  type(att_decomposer) :: att_mesh_decomposer_global
contains

  subroutine init_decomposer(this, nx, ny)
    class(att_decomposer), intent(inout) :: this
    integer, intent(in) :: nx, ny
    integer :: f1, f2, istart, iend
    integer, dimension(mysize,2) :: loc_ix, loc_iy
    integer :: max_rank_x, max_rank_y, myrank_value
    
    call close_factors(nx, ny, mysize, this%glob_px, this%glob_py)

    if (nx < this%glob_px .or. ny < this%glob_py) then
      call exit_MPI(myrank,'Error: Number of processors is larger than the number of cells')
    end if

    allocate(this%glob_ix(mysize, 2))
    allocate(this%glob_iy(mysize, 2))
    this%glob_ix = 0
    this%glob_iy = 0
    loc_ix = 0
    loc_iy = 0
    this%loc_ix_start = 1 + (mod(myrank, this%glob_px)) * (nx / this%glob_px)
    this%loc_ix_end = this%loc_ix_start + (nx / this%glob_px) - 1
    this%loc_iy_start = 1 + myrank / this%glob_px * (ny / this%glob_py)
    this%loc_iy_end = this%loc_iy_start + (ny / this%glob_py) - 1

    ! Adjust the end indices if am%n_xyz is odd
    if (mod(nx, this%glob_px) /= 0) then
      myrank_value = mod(myrank, this%glob_px)
      call max_all_all_i(myrank_value, max_rank_x)
      if (myrank_value==max_rank_x) this%loc_ix_end = nx
    end if
    if (mod(ny, this%glob_py) /= 0) then
      myrank_value = myrank / this%glob_px
      call max_all_all_i(myrank_value, max_rank_y)
      if (myrank_value==max_rank_y) this%loc_iy_end = ny
    end if      

    ! local number of cells
    loc_ix(myrank+1, 1) = this%loc_ix_start
    loc_ix(myrank+1, 2) = this%loc_ix_end
    loc_iy(myrank+1, 1) = this%loc_iy_start
    loc_iy(myrank+1, 2) = this%loc_iy_end
    this%loc_nx = this%loc_ix_end - this%loc_ix_start + 1
    this%loc_ny = this%loc_iy_end - this%loc_iy_start + 1
    call sum_all(loc_ix, this%glob_ix, mysize, 2)
    call sum_all(loc_iy, this%glob_iy, mysize, 2)

  end subroutine init_decomposer

  function partial_decompose_mesh(this, data) result(partial_data)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: data
    real(kind=dp), dimension(:,:,:), allocatable :: partial_data
    integer :: nz

    nz = size(data, 3)
    allocate(partial_data(this%loc_nx, this%loc_ny, nz))
    partial_data = data(this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :)

  end function partial_decompose_mesh

  subroutine collect_mesh(this, data, gathered_data)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: data
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_data
    real(kind=dp), dimension(:,:,:), intent(inout) :: gathered_data
    integer :: i, ierr, nx, ny, nz, tag=100

    nz = size(data, 3)
    if (myrank == 0) then
      gathered_data(this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :) = data
      do i = 1, mysize - 1
        nx = this%glob_ix(i+1, 2) - this%glob_ix(i+1, 1) + 1
        ny = this%glob_iy(i+1, 2) - this%glob_iy(i+1, 1) + 1
        allocate(tmp_data(nx, ny, nz))
        call recv_dp(tmp_data, nx * ny *nz, i, tag)
        gathered_data(this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2), :) = tmp_data
        deallocate(tmp_data)
      end do
    else
      allocate(tmp_data(this%loc_nx, this%loc_ny, nz))
      tmp_data = data
      call send_dp(tmp_data, this%loc_nx * this%loc_ny * nz, 0, tag)
      deallocate(tmp_data)
    end if
  end subroutine collect_mesh

  subroutine collect_grid(this, data, gathered_data)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: data
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_data
    real(kind=dp), dimension(:,:,:), intent(inout) :: gathered_data
    integer :: i, ierr, nx, ny, np, tag=100

    np = size(data, 1)
    if (myrank == 0) then
      gathered_data(:, this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end) = data
      do i = 1, mysize - 1
        nx = this%glob_ix(i+1, 2) - this%glob_ix(i+1, 1) + 1
        ny = this%glob_iy(i+1, 2) - this%glob_iy(i+1, 1) + 1
        allocate(tmp_data(np, nx, ny))
        call recv_dp(tmp_data, nx * ny * np, i, tag)
        gathered_data(:, this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2)) = tmp_data
        deallocate(tmp_data)
      end do
    else
      allocate(tmp_data(np, this%loc_nx, this%loc_ny))
      tmp_data = data
      call send_dp(tmp_data, this%loc_nx * this%loc_ny * np, 0, tag)
      deallocate(tmp_data)
    end if

  end subroutine collect_grid

  subroutine collect_sen(this, sen_vs, sen_vp, sen_rho, gathered_sen_vs, gathered_sen_vp, gathered_sen_rho)
    class(att_decomposer), intent(in) :: this
    real(kind=dp), dimension(:,:,:,:), intent(in) :: sen_vs, sen_vp, sen_rho
    real(kind=dp), dimension(:,:,:,:), allocatable :: tmp_sen_vs, tmp_sen_vp, tmp_sen_rho
    real(kind=dp), dimension(:,:,:,:), intent(inout) :: gathered_sen_vs, gathered_sen_vp, gathered_sen_rho
    integer :: i, ierr, nx, ny, nz, tag=100, np

    np = size(sen_vs, 1)
    nz = size(sen_vs, 4)
    if (myrank == 0) then
      gathered_sen_vs(:,this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :) = sen_vs
      gathered_sen_vp(:,this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :) = sen_vp
      gathered_sen_rho(:,this%loc_ix_start:this%loc_ix_end, this%loc_iy_start:this%loc_iy_end, :) = sen_rho
      do i = 1, mysize - 1
        nx = this%glob_ix(i+1, 2) - this%glob_ix(i+1, 1) + 1
        ny = this%glob_iy(i+1, 2) - this%glob_iy(i+1, 1) + 1
        allocate(tmp_sen_vs(np,nx, ny, nz))
        allocate(tmp_sen_vp(np,nx, ny, nz))
        allocate(tmp_sen_rho(np,nx, ny, nz))
        call recv_dp(tmp_sen_vs, np * nx * ny * nz, i, tag)
        call recv_dp(tmp_sen_vp, np * nx * ny * nz, i, tag)
        call recv_dp(tmp_sen_rho, np * nx * ny * nz, i, tag)
        gathered_sen_vs(:, this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2), :) = tmp_sen_vs
        gathered_sen_vp(:, this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2), :) = tmp_sen_vp
        gathered_sen_rho(:, this%glob_ix(i+1, 1):this%glob_ix(i+1, 2), this%glob_iy(i+1, 1):this%glob_iy(i+1, 2), :) = tmp_sen_rho
        deallocate(tmp_sen_vs)
        deallocate(tmp_sen_vp)
        deallocate(tmp_sen_rho)
      end do
    else
      allocate(tmp_sen_vs(np, this%loc_nx, this%loc_ny, nz))
      allocate(tmp_sen_vp(np, this%loc_nx, this%loc_ny, nz))
      allocate(tmp_sen_rho(np, this%loc_nx, this%loc_ny, nz))
      tmp_sen_vs = sen_vs
      tmp_sen_vp = sen_vp
      tmp_sen_rho = sen_rho
      call send_dp(tmp_sen_vs, np * this%loc_nx * this%loc_ny * nz, 0, tag)
      call send_dp(tmp_sen_vp, np * this%loc_nx * this%loc_ny * nz, 0, tag)
      call send_dp(tmp_sen_rho, np * this%loc_nx * this%loc_ny * nz, 0, tag)
      deallocate(tmp_sen_vs)
      deallocate(tmp_sen_vp)
      deallocate(tmp_sen_rho)
    end if  

  end subroutine collect_sen

  subroutine close_factors(nx, ny, num, f1o, f2o)
    integer, intent(in) :: num, nx, ny
    integer, intent(out) :: f1o, f2o
    integer :: f1, f2, i
    real :: dif0,dd,dif1
    
    dif0 = real(nx)/real(ny)
    f1 = 1; f2 = num
    f1o = f1; f2o = f2
    dd = abs(dif0 - real(f1)/real(f2))
    do i = 2, num
      if (mod(num, i) == 0) then
        f1 = i
        f2 = num / i
        dif1 = real(f1)/real(f2)
        if (abs(dif0-dif1) < dd) then
          dd = abs(dif0-dif1)
          f1o = f1; f2o=f2
          cycle
        else
          return
        endif
      end if
    end do
  end subroutine close_factors
end module