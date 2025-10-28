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
!                       May 2024, Add L-BFGS and CG optimization
!                       Jun 2024, Add Hager-Zhang formula for CG
!                       Jan 2025, change dims for Alpha_beta_rho
!
!=====================================================================
module optimize
  use hdf5_interface
  use shared_par
  use para, ap => att_para_global
  use model, am => att_model_global
  use utils

  implicit none
contains
  subroutine get_lbfgs_direction(iter, direction)
    integer, intent(in) :: iter
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: direction
    real(kind=dp), dimension(:,:,:,:), allocatable :: gradient0,gradient1,model0,model1,&
                                                    q_vector,r_vector
    real(kind=dp), dimension(:,:,:,:,:), allocatable :: gradient_diff, model_diff
    real(kind=dp), dimension(:), allocatable :: p, a
    real(kind=dp) :: p_k_up_sum, p_k_down_sum, p_k, b
    integer :: istore, i, dims(4), iter_store, nstore
    integer, dimension(:), allocatable :: idx_iter

    iter_store = iter-m_store
    if (iter_store <= iter_start) iter_store = iter_start
    nstore = iter-iter_store

    call get_gradient(iter, q_vector)
    dims = shape(q_vector)
    allocate(gradient_diff(nstore, nker, dims(2), dims(3), dims(4)))
    allocate(model_diff(nstore, nker, dims(2), dims(3), dims(4)))
    idx_iter = zeros(nstore)
    p = zeros(nstore)
    a = zeros(nstore)
    i = 0
    do istore = iter-1,iter_store,-1
      i = i+1
      idx_iter(i) = istore
      call get_gradient(istore, gradient0)
      call get_gradient(istore+1, gradient1)
      call get_model(istore, model0)
      call get_model(istore+1, model1)

      model_diff(i,:,:,:,:) = model1 - model0
      gradient_diff(i,:,:,:,:) = gradient1 - gradient0

      p(i) = 1/sum(model_diff(i,:,:,:,:)*gradient_diff(i,:,:,:,:))
      a(i) = sum(model_diff(i,:,:,:,:)*q_vector)*p(i)
      q_vector = q_vector - a(i)*gradient_diff(i,:,:,:,:)
    enddo
    p_k_up_sum = sum(model_diff(1,:,:,:,:)*gradient_diff(1,:,:,:,:))
    p_k_down_sum = sum(gradient_diff(1,:,:,:,:)*gradient_diff(1,:,:,:,:))
    p_k = p_k_up_sum/p_k_down_sum
    r_vector = p_k*q_vector

    do istore = iter_store,iter-1
      i = find_loc(idx_iter, istore)
      b = sum(gradient_diff(i,:,:,:,:)*r_vector)*p(i)
      r_vector = r_vector + model_diff(i,:,:,:,:)*(a(i)-b)
    enddo
    direction = -1.0_dp * r_vector

  end subroutine get_lbfgs_direction

  subroutine get_cg_direction(iter, direction)
    ! get the conjugate gradient direction using Hager-Zhang formula
    integer, intent(in) :: iter
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: direction
    real(kind=dp), dimension(:,:,:,:), allocatable :: gradient0,gradient1,direction0,grad_diff
    real(kind=dp) :: beta

    call get_gradient(iter-1, gradient0)
    call get_gradient(iter, gradient1)
    call get_direction(iter-1, direction0)
    grad_diff = gradient1 - gradient0
    beta = sum(grad_diff * grad_diff) / sum(grad_diff * direction0)
    beta = sum(gradient1*(grad_diff-2*direction0*beta)) / sum(direction0*grad_diff)
    ! beta = sum(gradient1*(gradient1-gradient0))/sum(gradient0*gradient0)
    direction = -gradient1 + beta*direction0
    
  end subroutine get_cg_direction

  subroutine get_model(iter, model)
    integer, intent(in) :: iter
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: model
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_model
    character(len=MAX_NAME_LEN) :: key_name

    write(key_name, '("/vs_",I3.3)') iter
    ! call h5read(model_fname, key_name, tmp_model)
    call am%modfp%get(key_name, tmp_model)
    tmp_model = transpose_3(tmp_model)
    model = zeros( nker,size(tmp_model, 1), size(tmp_model, 2), size(tmp_model, 3))
    model(1,:,:,:) = tmp_model
    if (ap%inversion%use_alpha_beta_rho) then
      ! read vp model
      deallocate(tmp_model)
      write(key_name, '("/vp_",I3.3)') iter
      ! call h5read(model_fname, key_name, tmp_model)
      call am%modfp%get(key_name, tmp_model)
      model(2,:,:,:) = transpose_3(tmp_model)

      ! read rho model
      deallocate(tmp_model)
      write(key_name, '("/rho_",I3.3)') iter
      ! call h5read(model_fname, key_name, tmp_model)
      call am%modfp%get(key_name, tmp_model)
      model(3,:,:,:) = transpose_3(tmp_model)
    endif
    model = log(model)

  end subroutine get_model

  subroutine get_gradient(iter, gradient)
    integer, intent(in) :: iter
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: gradient
    character(len=MAX_NAME_LEN) :: key_name

    write(key_name, '("/gradient_",I3.3)') iter
    ! call h5read(model_fname, key_name, gradient)
    call am%modfp%get(key_name, gradient)
    gradient = transpose_4(gradient)

  end subroutine get_gradient

  subroutine get_direction(iter, direction)
    integer, intent(in) :: iter
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: direction
    character(len=MAX_NAME_LEN) :: key_name

    write(key_name, '("/direction_",I3.3)') iter
    ! call h5read(model_fname, key_name, direction)
    call am%modfp%get(key_name, direction)
    direction = transpose_4(direction)

  end subroutine
end module optimize

