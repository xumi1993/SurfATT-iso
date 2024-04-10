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
! This file include the subroutines for solving the eikonal equation 
!     and its adjoint equation in 2D domain. These codes are cridited
!     by Chen Jing at NTU.
! This file also include subroutines for conversion between muti-grid
!      and forward mesh for optimization, which credits to Ping Tong at NTU.
!
!=====================================================================

! FILE: eikon_solver_2d.f90
subroutine FSM_UW_PS_cart_2d(xx,yy,nx,ny,a,b,c,T,fun,x0,y0,u)
    integer :: nx,ny
    double precision :: xx(nx),yy(ny),a(nx,ny),b(nx,ny),c(nx,ny),T(nx,ny),fun(nx,ny),x0,y0,ischange(nx,ny),u(nx,ny)
    double precision :: dx,dy,T0(nx,ny),T0x(nx,ny),T0y(nx,ny),a0,b0,c0,fun0
    double precision :: tau(nx,ny)
    integer :: iix,iiy,iter,i_sweep, x_id1,x_id2,y_id1,y_id2,xdirec,ydirec
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nx,ny)
    integer,parameter :: MaxIter = 2000
    double precision,parameter :: tol=(10.0)**(-4),eps=(10.0)**(-14)

    integer :: i_case,i_cand,i_solution,count_cand
    double precision :: px,qx,py,qy
    double precision :: eqn_a,eqn_b,eqn_c,Delta,tmp_tau,dis
    double precision :: characteristic_x, characteristic_y, Tx, Ty
    logical :: is_causility
    double precision :: canditate(50)

!f2py   intent(in) :: xx, yy
!f2py   intent(hide), depend(xx) :: nx = shape(xx, 0)
!f2py   intent(hide), depend(yy) :: ny = shape(yy, 0)
!f2py   intent(in) :: a, b, c, u, fun
!f2py   intent(in) :: x0, y0
!f2py   intent(out) :: T

    dx=xx(2)-xx(1); dy=yy(2)-yy(1);


    idx0=floor((x0-xx(1))/dx+1); idy0=floor((y0-yy(1))/dy+1); 
    r1 = min(1.0,(x0-xx(idx0))/dx); r2 = min(1.0,(y0-yy(idy0))/dy); 

    a0=(1-r1)*(1-r2)*a(idx0,idy0)+(1-r1)*r2*a(idx0,idy0+1) &
    & +r1*(1-r2)*a(idx0+1,idy0)+r1*r2*a(idx0+1,idy0+1) 

    b0=(1-r1)*(1-r2)*b(idx0,idy0)+(1-r1)*r2*b(idx0,idy0+1) &
    & +r1*(1-r2)*b(idx0+1,idy0)+r1*r2*b(idx0+1,idy0+1) 

    c0=(1-r1)*(1-r2)*c(idx0,idy0)+(1-r1)*r2*c(idx0,idy0+1) &
    & +r1*(1-r2)*c(idx0+1,idy0)+r1*r2*c(idx0+1,idy0+1) 

    fun0=(1-r1)*(1-r2)*fun(idx0,idy0)+(1-r1)*r2*fun(idx0,idy0+1) &
    & +r1*(1-r2)*fun(idx0+1,idy0)+r1*r2*fun(idx0+1,idy0+1) 


    ! solve T0 = sqrt((b0(x-x0)^2+a0(y-y0)^2+2c0(x-x0)(y-y0))/(a0b0-c0^2))
    ! write(*,*) fun0,b0/(a0*b0-c0**2),a0/(a0*b0-c0**2),2*c0/(a0*b0-c0**2)
    do iix=1,nx
        do iiy=1,ny
            T0(iix,iiy) = fun0*sqrt((b0/(a0*b0-c0**2)*(xx(iix)-x0)**2 + a0/(a0*b0-c0**2)*(yy(iiy)-y0)**2 &
                        & + 2*c0/(a0*b0-c0**2)*(xx(iix)-x0)*(yy(iiy)-y0)))
            if (T0(iix,iiy) .eq. 0) then
                T0x(iix,iiy) = 0
                T0y(iix,iiy) = 0
            else
                T0x(iix,iiy) = fun0**2*(b0/(a0*b0-c0**2)*(xx(iix)-x0)+c0/(a0*b0-c0**2)*(yy(iiy)-y0))/T0(iix,iiy)
                T0y(iix,iiy) = fun0**2*(a0/(a0*b0-c0**2)*(yy(iiy)-y0)+c0/(a0*b0-c0**2)*(xx(iix)-x0))/T0(iix,iiy)
            end if

            if ( abs((xx(iix)-x0)/dx)<=1 .and. abs((yy(iiy)-y0)/dy)<=1) then
                tau(iix,iiy) = 1 
                ischange(iix,iiy)=0
                if (iix==1 .or. iix==nx .or. iiy==1 .or. iiy==ny) then
                    write(*,*) 'source on the boundary, mesh error'
                    stop
                end if
            else
                tau(iix,iiy) = 10
                ischange(iix,iiy)=1
            end if


        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
    ! do iter =1,1
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do i_sweep = 1,4
        ! do i_sweep = 1,1
            select case(i_sweep)
                case (1)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (2)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1
                case (3)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (4)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1    
            end select


            ! iter 1 x: 1 -> Nx, y: 1 -> Ny
            ! iter 2 x: 1 -> Nx, y: Ny -> 1 
            ! iter 3 x: Nx -> 1, y: 1 -> Ny 
            ! iter 4 x: Nx -> 1, y: Ny -> 1 
            do iix=x_id1,x_id2,xdirec
                do iiy=y_id1,y_id2,ydirec
                    if(ischange(iix,iiy)==1) then
    !--------------------------------------- calculate stencil start
                        count_cand = 1
                        do i_case=1,4       
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)    !  -1, -1
                                if (iix == 1 .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)     
                            case (2)    ! -1, +1
                                if (iix == 1 .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1)   
                            case (3)    !  +1, -1
                                if (iix == Nx .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx = +T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                            case (4)    ! 右上  +1, +1
                                if (iix == Nx .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                            end select
                            
                            eqn_a = a(iix,iiy)*px**2 + b(iix,iiy)*py**2 -2*c(iix,iiy)*px*py
                            eqn_b = 2*a(iix,iiy)*px*qx + 2*b(iix,iiy)*py*qy - 2*c(iix,iiy)*(px*qy+py*qx)
                            eqn_c = a(iix,iiy)*qx**2 + b(iix,iiy)*qy**2 -2*c(iix,iiy)*qx*qy-fun(iix,iiy)**2
                            Delta = eqn_b**2-4*eqn_a*eqn_c
                            

                            if (Delta>=0) then  !   one or two solutions
                                ! check the causality condition: the characteristic passing through (iix,iiy) is in between used two sides
                                do i_solution=1,2
                                    select case (i_solution)
                                    case(1)
                                        tmp_tau = (-eqn_b + sqrt(Delta))/(2*eqn_a)
                                    case(2)
                                        tmp_tau = (-eqn_b - sqrt(Delta))/(2*eqn_a)
                                    end select

                                    ! characteristic is
                                    Tx = px*tmp_tau+qx
                                    Ty = py*tmp_tau+qy
                                    characteristic_x = a(iix,iiy)*Tx - c(iix,iiy)*Ty 
                                    characteristic_y = b(iix,iiy)*Ty - c(iix,iiy)*Tx 

                                    is_causility = .false.
                                    select case (i_case)
                                    case (1)    
                                        if (characteristic_x >= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (2)   
                                        if (characteristic_x >= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    case (3)  
                                        if (characteristic_x <= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (4)   
                                        if (characteristic_x <= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    end select

                                    ! if satisfying the causility condition, retain it as a canditate solution
                                    if(is_causility) then
                                        canditate(count_cand) = tmp_tau
                                        count_cand = count_cand + 1
                                    end if


                                end do
                            end if

                        end do

                        do i_case=1,4   
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; 
                                ! (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)    
                                if (iix == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy) &
                                & .and. tmp_tau > tau(iix-1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy) &
                                & .and. tmp_tau > tau(iix-1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (2)   
                                if (iix == Nx) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy) &
                                & .and. tmp_tau > tau(iix+1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy) &
                                & .and. tmp_tau > tau(iix+1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if

                            case (3)    
                                if (iiy == 1) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1) &
                                & .and. tmp_tau > tau(iix,iiy-1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1) &
                                & .and. tmp_tau > tau(iix,iiy-1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (4)    
                                if (iiy == Ny) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1) &
                                & .and. tmp_tau > tau(iix,iiy+1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1) &
                                & .and. tmp_tau > tau(iix,iiy+1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            end select

                        end do

                        do i_cand=1,count_cand-1
                            tau(iix,iiy) = min(tau(iix,iiy),canditate(i_cand))
                        end do
    !--------------------------------------- calculate stencil end


                    end if
                end do
            end do
                
                
                
        end do

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy)
                Linf_dif=max(Linf_dif,abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy))               
            end do
        end do

        do iix=1,nx
            do iiy=1,ny
                L1_err=L1_err+abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy))
                Linf_err=max(Linf_err,abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy)))            
            end do
        end do
        L1_dif = L1_dif/(nx)/(ny)
        L1_err = L1_err/(nx)/(ny)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            ! write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
            ! write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if
        ! write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        ! write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau*T0

end subroutine

subroutine FSM_WENO3_PS_2d(xx,yy,nx,ny,a,b,c,T,fun,x0,y0,u)
    integer :: nx,ny
    double precision :: xx(nx),yy(ny),a(nx,ny),b(nx,ny),c(nx,ny),T(nx,ny),fun(nx,ny),x0,y0,ischange(nx,ny),u(nx,ny)
    double precision :: dx,dy,T0(nx,ny),T0x(nx,ny),T0y(nx,ny),a0,b0,c0,fun0
    double precision :: tau(nx,ny),px1,px2,py1,py2,tpT,Htau,wx1,wx2,wy1,wy2
    integer :: iix,iiy,iter,xdirec,ydirec
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nx,ny)
    integer,parameter :: MaxIter = 1000
    double precision,parameter :: tol=(10.0)**(-6),eps=(10.0)**(-14)

!f2py   intent(in) :: xx, yy
!f2py   intent(hide) :: nx = shape(xx, 0)
!f2py   intent(hide) :: ny = shape(yy, 0)
!f2py   intent(in) :: a, b, c, u, fun
!f2py   intent(in) :: x0, y0
!f2py   intent(out) :: T

    ! ------------------------ create mesh ------------------------ 
    dx=xx(2)-xx(1); dy=yy(2)-yy(1);

    ! ------------------------ create T0 ------------------------ 

    ! Parameter discretization at the seismic source.
    idx0=floor((x0-xx(1))/dx+1); idy0=floor((y0-yy(1))/dy+1); 
    r1 = min(1.0,(x0-xx(idx0))/dx); r2 = min(1.0,(y0-yy(idy0))/dy); 

    a0=(1-r1)*(1-r2)*a(idx0,idy0)+(1-r1)*r2*a(idx0,idy0+1) &
    & +r1*(1-r2)*a(idx0+1,idy0)+r1*r2*a(idx0+1,idy0+1) 

    b0=(1-r1)*(1-r2)*b(idx0,idy0)+(1-r1)*r2*b(idx0,idy0+1) &
    & +r1*(1-r2)*b(idx0+1,idy0)+r1*r2*b(idx0+1,idy0+1) 

    c0=(1-r1)*(1-r2)*c(idx0,idy0)+(1-r1)*r2*c(idx0,idy0+1) &
    & +r1*(1-r2)*c(idx0+1,idy0)+r1*r2*c(idx0+1,idy0+1) 

    fun0=(1-r1)*(1-r2)*fun(idx0,idy0)+(1-r1)*r2*fun(idx0,idy0+1) &
    & +r1*(1-r2)*fun(idx0+1,idy0)+r1*r2*fun(idx0+1,idy0+1) 


    ! solve T0 = sqrt((b0(x-x0)^2+a0(y-y0)^2+2c0(x-x0)(y-y0))/(a0b0-c0^2))
    !write(*,*) fun0,b0/(a0*b0-c0**2),a0/(a0*b0-c0**2),2*c0/(a0*b0-c0**2)
    do iix=1,nx
        do iiy=1,ny
            T0(iix,iiy) = fun0*sqrt((b0/(a0*b0-c0**2)*(xx(iix)-x0)**2 + a0/(a0*b0-c0**2)*(yy(iiy)-y0)**2 &
                        & + 2*c0/(a0*b0-c0**2)*(xx(iix)-x0)*(yy(iiy)-y0)))
            if (T0(iix,iiy) .eq. 0) then
                T0x(iix,iiy) = 0
                T0y(iix,iiy) = 0
            else
                T0x(iix,iiy) = fun0**2*(b0/(a0*b0-c0**2)*(xx(iix)-x0)+c0/(a0*b0-c0**2)*(yy(iiy)-y0))/T0(iix,iiy)
                T0y(iix,iiy) = fun0**2*(a0/(a0*b0-c0**2)*(yy(iiy)-y0)+c0/(a0*b0-c0**2)*(xx(iix)-x0))/T0(iix,iiy)
            end if

            if ( abs((xx(iix)-x0)/dx)<=2 .and. abs((yy(iiy)-y0)/dy)<=2) then
                tau(iix,iiy) = 0  
                ischange(iix,iiy)=0
                if (iix==1 .or. iix==nx .or. iiy==1 .or. iiy==ny) then
                    write(*,*) 'source on the boundary, mesh error'
                    stop
                end if
            else
                tau(iix,iiy) = 20
                ischange(iix,iiy)=1
            end if

        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do xdirec = -1,1,2
            do ydirec = -1,1,2
                ! iter 1 x: 1 -> Nx, y: 1 -> Ny
                ! iter 2 x: 1 -> Nx, y: Ny -> 1 
                ! iter 3 x: Nx -> 1, y: 1 -> Ny 
                ! iter 4 x: Nx -> 1, y: Ny -> 1 
                do iix=nint(0.5+Nx/2.0+(Nx/2.0-1.5)*xdirec),nint(0.5+Nx/2.0+(-Nx/2.0+1.5)*xdirec),-xdirec
                    do iiy=nint(0.5+Ny/2.0+(Ny/2.0-1.5)*ydirec),nint(0.5+Ny/2.0+(-Ny/2.0+1.5)*ydirec),-ydirec
                        if(ischange(iix,iiy)==1) then
                            
                            sigx=sqrt(a(iix,iiy));sigy=sqrt(b(iix,iiy))
                            coe=1.0/((sigx/dx)+(sigy/dy)) 

                            if(iix==2) then
                                px1=(tau(iix,iiy)-tau(iix-1,iiy))/dx; 
                                wx2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix+1,iiy)+tau(iix+2,iiy))**2)/ &
                                            & (eps+(tau(iix-1,iiy)-2*tau(iix,iiy)+tau(iix+1,iiy))**2))**2)
                                px2=(1-wx2)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx2*(-3*tau(iix,iiy)+4*tau(iix+1,iiy)-tau(iix+2,iiy))/2/dx;
                            elseif (iix==Nx-1) then
                                wx1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix-1,iiy)+tau(iix-2,iiy))**2)/ &
                                            & (eps+(tau(iix+1,iiy)-2*tau(iix,iiy)+tau(iix-1,iiy))**2))**2)
                                px1=(1-wx1)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx1*(3*tau(iix,iiy)-4*tau(iix-1,iiy)+tau(iix-2,iiy))/2/dx;
                                px2=(tau(iix+1,iiy)-tau(iix,iiy))/dx; 
                            else
                                wx1=1.0/(1.0+2*((eps+(tau(iix,iiy)-2*tau(iix-1,iiy)+tau(iix-2,iiy))**2)/ &
                                            & (eps+(tau(iix+1,iiy)-2*tau(iix,iiy)+tau(iix-1,iiy))**2))**2)
                                px1=(1.0-wx1)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx1*(3*tau(iix,iiy)-4*tau(iix-1,iiy)+tau(iix-2,iiy))/2/dx;
                                wx2=1.0/(1.0+2*((eps+(tau(iix,iiy)-2*tau(iix+1,iiy)+tau(iix+2,iiy))**2)/ &
                                            & (eps+(tau(iix-1,iiy)-2*tau(iix,iiy)+tau(iix+1,iiy))**2))**2)
                                px2=(1.0-wx2)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx2*(-3*tau(iix,iiy)+4*tau(iix+1,iiy)-tau(iix+2,iiy))/2/dx;
                            end if


                            if(iiy==2) then
                                py1=(tau(iix,iiy)-tau(iix,iiy-1))/dy; 
                                wy2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy+1)+tau(iix,iiy+2))**2)/ &
                                            & (eps+(tau(iix,iiy-1)-2*tau(iix,iiy)+tau(iix,iiy+1))**2))**2)
                                py2=(1-wy2)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy2*(-3*tau(iix,iiy)+4*tau(iix,iiy+1)-tau(iix,iiy+2))/2/dy;
                            elseif (iiy==Ny-1) then
                                wy1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy-1)+tau(iix,iiy-2))**2)/ &
                                            & (eps+(tau(iix,iiy+1)-2*tau(iix,iiy)+tau(iix,iiy-1))**2))**2)
                                py1=(1-wy1)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy1*(3*tau(iix,iiy)-4*tau(iix,iiy-1)+tau(iix,iiy-2))/2/dy;
                                py2=(tau(iix,iiy+1)-tau(iix,iiy))/dy;
                            else
                                wy1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy-1)+tau(iix,iiy-2))**2)/ &
                                            & (eps+(tau(iix,iiy+1)-2*tau(iix,iiy)+tau(iix,iiy-1))**2))**2)
                                py1=(1-wy1)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy1*(3*tau(iix,iiy)-4*tau(iix,iiy-1)+tau(iix,iiy-2))/2/dy;
                                wy2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy+1)+tau(iix,iiy+2))**2)/ &
                                            & (eps+(tau(iix,iiy-1)-2*tau(iix,iiy)+tau(iix,iiy+1))**2))**2)
                                py2=(1-wy2)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy2*(-3*tau(iix,iiy)+4*tau(iix,iiy+1)-tau(iix,iiy+2))/2/dy;                            
                            end if
                
                            Htau=sqrt( a(iix,iiy)*((px1+px2)/2)**2 + b(iix,iiy)*((py1+py2)/2)**2 &
                                    & -2*c(iix,iiy)*(px1+px2)*(py1+py2)/4 &
                                    & +2*(a(iix,iiy)*T0x(iix,iiy)-c(iix,iiy)*T0y(iix,iiy))*(px1+px2)/2 &
                                    & +2*(b(iix,iiy)*T0y(iix,iiy)-c(iix,iiy)*T0x(iix,iiy))*(py1+py2)/2 &
                                    & +a(iix,iiy)*T0x(iix,iiy)**2+b(iix,iiy)*T0y(iix,iiy)**2 &
                                    & -2*c(iix,iiy)*T0x(iix,iiy)*T0y(iix,iiy) )

                            !sigx=min(1.0,abs(a(iix,iiy)*((px1+px2)/2+T0x(iix,iiy))-c(iix,iiy)*(T0y(iix,iiy)+(py1+py2)/2))/Htau*2)
                            !sigy=min(1.0,abs(b(iix,iiy)*((py1+py2)/2+T0y(iix,iiy))-c(iix,iiy)*(T0x(iix,iiy)+(px1+px2)/2))/Htau*2)
                            
                            !sigx=1!abs(a(iix,iiy)*((px1+px2)/2+T0x(iix,iiy))-c(iix,iiy)*(T0y(iix,iiy)+(py1+py2)/2))/Htau*2
                            !sigy=1!abs(b(iix,iiy)*((py1+py2)/2+T0y(iix,iiy))-c(iix,iiy)*(T0x(iix,iiy)+(px1+px2)/2))/Htau*2
                            
                            tpT=coe*(fun(iix,iiy)-Htau) + coe*(sigx*(px2-px1)/2+sigy*(py2-py1)/2)+tau(iix,iiy);
                                                        
                            tau(iix,iiy)=tpT

                            

                        end if

                    end do
                end do
                
                do iiy=1,ny
                    tau(1,iiy) = max(2*tau(2,iiy)-tau(3,iiy),tau(3,iiy))
                    tau(nx,iiy) = max(2*tau(nx-1,iiy)-tau(nx-2,iiy),tau(nx-2,iiy))
                end do
                do iix=1,nx
                    tau(iix,1) = max(2*tau(iix,2)-tau(iix,3),tau(iix,3))
                    tau(iix,ny) = max(2*tau(iix,ny-1)-tau(iix,ny-2),tau(iix,ny-2))
                end do
                
            end do
        end do

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(tau(iix,iiy)-tau_old(iix,iiy))*dx*dy
                Linf_dif=max(Linf_dif,abs(tau(iix,iiy)-tau_old(iix,iiy)))               
            end do
        end do

        do iix=5,nx-4
            do iiy=5,ny-4
                L1_err=L1_err+abs(tau(iix,iiy)+T0(iix,iiy)-u(iix,iiy))
                Linf_err=max(Linf_err,abs(tau(iix,iiy)+T0(iix,iiy)-u(iix,iiy)))            
            end do
        end do

        L1_err = L1_err/(nx-4)/(ny-4)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
        !    write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if
        !write(*,'(a,f10.7,a,f10.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        !write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau+T0

end subroutine

! 2D Eikonal, more stable
! fun is slowness field, u is not important
subroutine FSM_WENO3_PS_multi_2d(xx,yy,nx,ny,a,b,c,T,fun,x0,y0,u)
    integer :: nx,ny
    double precision :: xx(nx),yy(ny),a(nx,ny),b(nx,ny),c(nx,ny),T(nx,ny),fun(nx,ny),x0,y0,ischange(nx,ny),u(nx,ny)
    double precision :: dx,dy,T0(nx,ny),T0x(nx,ny),T0y(nx,ny),a0,b0,c0,fun0
    double precision :: tau(nx,ny),px1,px2,py1,py2,tpT,Htau,wx1,wx2,wy1,wy2
    integer :: iix,iiy,iter,xdirec,ydirec
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nx,ny)
    integer,parameter :: MaxIter = 1000
    double precision,parameter :: tol=(10.0)**(-6),eps=(10.0)**(-14)

!f2py   intent(in) :: xx, yy
!f2py   intent(hide), depend(xx) :: nx = shape(xx, 0)
!f2py   intent(hide), depend(yy) :: ny = shape(yy, 0)
!f2py   intent(in) :: a, b, c, u, fun
!f2py   intent(in) :: x0, y0
!f2py   intent(out) :: T
    ! ------------------------ create mesh ------------------------ 
    dx=xx(2)-xx(1); dy=yy(2)-yy(1);

    ! ------------------------ create T0 ------------------------ 

    ! Parameter discretization at the seismic source.
    idx0=floor((x0-xx(1))/dx+1); idy0=floor((y0-yy(1))/dy+1); 
    r1 = min(1.0,(x0-xx(idx0))/dx); r2 = min(1.0,(y0-yy(idy0))/dy); 

    a0=(1-r1)*(1-r2)*a(idx0,idy0)+(1-r1)*r2*a(idx0,idy0+1) &
    & +r1*(1-r2)*a(idx0+1,idy0)+r1*r2*a(idx0+1,idy0+1) 

    b0=(1-r1)*(1-r2)*b(idx0,idy0)+(1-r1)*r2*b(idx0,idy0+1) &
    & +r1*(1-r2)*b(idx0+1,idy0)+r1*r2*b(idx0+1,idy0+1) 

    c0=(1-r1)*(1-r2)*c(idx0,idy0)+(1-r1)*r2*c(idx0,idy0+1) &
    & +r1*(1-r2)*c(idx0+1,idy0)+r1*r2*c(idx0+1,idy0+1) 

    fun0=(1-r1)*(1-r2)*fun(idx0,idy0)+(1-r1)*r2*fun(idx0,idy0+1) &
    & +r1*(1-r2)*fun(idx0+1,idy0)+r1*r2*fun(idx0+1,idy0+1) 


    ! solve T0 = sqrt((b0(x-x0)^2+a0(y-y0)^2+2c0(x-x0)(y-y0))/(a0b0-c0^2))
    !write(*,*) fun0,b0/(a0*b0-c0**2),a0/(a0*b0-c0**2),2*c0/(a0*b0-c0**2)
    do iix=1,nx
        do iiy=1,ny
            T0(iix,iiy) = fun0*sqrt((b0/(a0*b0-c0**2)*(xx(iix)-x0)**2 + a0/(a0*b0-c0**2)*(yy(iiy)-y0)**2 &
                        & + 2*c0/(a0*b0-c0**2)*(xx(iix)-x0)*(yy(iiy)-y0)))
            if (T0(iix,iiy) .eq. 0) then
                T0x(iix,iiy) = 0
                T0y(iix,iiy) = 0
            else
                T0x(iix,iiy) = fun0**2*(b0/(a0*b0-c0**2)*(xx(iix)-x0)+c0/(a0*b0-c0**2)*(yy(iiy)-y0))/T0(iix,iiy)
                T0y(iix,iiy) = fun0**2*(a0/(a0*b0-c0**2)*(yy(iiy)-y0)+c0/(a0*b0-c0**2)*(xx(iix)-x0))/T0(iix,iiy)
            end if

            if ( abs((xx(iix)-x0)/dx)<=2 .and. abs((yy(iiy)-y0)/dy)<=2) then
                !Several points around the seismic source are directly assumed to have a constant velocity structure, for which an analytical solution is given, i.e., equal to T0.
                tau(iix,iiy) = 1  
                ischange(iix,iiy)=0
                if (iix==1 .or. iix==nx .or. iiy==1 .or. iiy==ny) then
                    write(*,*) 'source on the boundary, mesh error'
                    write(*,*) x0, y0
                    stop
                end if
            else
                tau(iix,iiy) = 1
                ischange(iix,iiy)=1
            end if
        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do xdirec = -1,1,2
            do ydirec = -1,1,2
                ! iter 1 x: 1 -> Nx, y: 1 -> Ny
                ! iter 2 x: 1 -> Nx, y: Ny -> 1 
                ! iter 3 x: Nx -> 1, y: 1 -> Ny 
                ! iter 4 x: Nx -> 1, y: Ny -> 1 
                do iix=nint(0.5+Nx/2.0+(Nx/2.0-1.5)*xdirec),nint(0.5+Nx/2.0+(-Nx/2.0+1.5)*xdirec),-xdirec
                    do iiy=nint(0.5+Ny/2.0+(Ny/2.0-1.5)*ydirec),nint(0.5+Ny/2.0+(-Ny/2.0+1.5)*ydirec),-ydirec
                        if(ischange(iix,iiy)==1) then
                            
                            sigx=sqrt(a(iix,iiy))*T0(iix,iiy);sigy=sqrt(b(iix,iiy))*T0(iix,iiy)
                            coe=1.0/((sigx/dx)+(sigy/dy)) 

                            if(iix==2) then
                                px1=(tau(iix,iiy)-tau(iix-1,iiy))/dx; 
                                wx2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix+1,iiy)+tau(iix+2,iiy))**2)/ &
                                            & (eps+(tau(iix-1,iiy)-2*tau(iix,iiy)+tau(iix+1,iiy))**2))**2)
                                px2=(1-wx2)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx2*(-3*tau(iix,iiy)+4*tau(iix+1,iiy)-tau(iix+2,iiy))/2/dx;
                            elseif (iix==Nx-1) then
                                wx1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix-1,iiy)+tau(iix-2,iiy))**2)/ &
                                            & (eps+(tau(iix+1,iiy)-2*tau(iix,iiy)+tau(iix-1,iiy))**2))**2)
                                px1=(1-wx1)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx1*(3*tau(iix,iiy)-4*tau(iix-1,iiy)+tau(iix-2,iiy))/2/dx;
                                px2=(tau(iix+1,iiy)-tau(iix,iiy))/dx; 
                            else
                                wx1=1.0/(1.0+2*((eps+(tau(iix,iiy)-2*tau(iix-1,iiy)+tau(iix-2,iiy))**2)/ &
                                            & (eps+(tau(iix+1,iiy)-2*tau(iix,iiy)+tau(iix-1,iiy))**2))**2)
                                px1=(1.0-wx1)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx1*(3*tau(iix,iiy)-4*tau(iix-1,iiy)+tau(iix-2,iiy))/2/dx;
                                wx2=1.0/(1.0+2*((eps+(tau(iix,iiy)-2*tau(iix+1,iiy)+tau(iix+2,iiy))**2)/ &
                                            & (eps+(tau(iix-1,iiy)-2*tau(iix,iiy)+tau(iix+1,iiy))**2))**2)
                                px2=(1.0-wx2)*(tau(iix+1,iiy)-tau(iix-1,iiy))/2/dx+ &
                                      & wx2*(-3*tau(iix,iiy)+4*tau(iix+1,iiy)-tau(iix+2,iiy))/2/dx;
                            end if


                            if(iiy==2) then
                                py1=(tau(iix,iiy)-tau(iix,iiy-1))/dy; 
                                wy2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy+1)+tau(iix,iiy+2))**2)/ &
                                            & (eps+(tau(iix,iiy-1)-2*tau(iix,iiy)+tau(iix,iiy+1))**2))**2)
                                py2=(1-wy2)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy2*(-3*tau(iix,iiy)+4*tau(iix,iiy+1)-tau(iix,iiy+2))/2/dy;
                            elseif (iiy==Ny-1) then
                                wy1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy-1)+tau(iix,iiy-2))**2)/ &
                                            & (eps+(tau(iix,iiy+1)-2*tau(iix,iiy)+tau(iix,iiy-1))**2))**2)
                                py1=(1-wy1)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy1*(3*tau(iix,iiy)-4*tau(iix,iiy-1)+tau(iix,iiy-2))/2/dy;
                                py2=(tau(iix,iiy+1)-tau(iix,iiy))/dy;
                            else
                                wy1=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy-1)+tau(iix,iiy-2))**2)/ &
                                            & (eps+(tau(iix,iiy+1)-2*tau(iix,iiy)+tau(iix,iiy-1))**2))**2)
                                py1=(1-wy1)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy1*(3*tau(iix,iiy)-4*tau(iix,iiy-1)+tau(iix,iiy-2))/2/dy;
                                wy2=1.0/(1+2*((eps+(tau(iix,iiy)-2*tau(iix,iiy+1)+tau(iix,iiy+2))**2)/ &
                                            & (eps+(tau(iix,iiy-1)-2*tau(iix,iiy)+tau(iix,iiy+1))**2))**2)
                                py2=(1-wy2)*(tau(iix,iiy+1)-tau(iix,iiy-1))/2/dy+ &
                                    & wy2*(-3*tau(iix,iiy)+4*tau(iix,iiy+1)-tau(iix,iiy+2))/2/dy;                            
                            end if
                
                            !Htau=sqrt( a(iix,iiy)*((px1+px2)/2)**2 + b(iix,iiy)*((py1+py2)/2)**2 &
                            !        & -2*c(iix,iiy)*(px1+px2)*(py1+py2)/4 &
                            !        & +2*(a(iix,iiy)*T0x(iix,iiy)-c(iix,iiy)*T0y(iix,iiy))*(px1+px2)/2 &
                            !        & +2*(b(iix,iiy)*T0y(iix,iiy)-c(iix,iiy)*T0x(iix,iiy))*(py1+py2)/2 &
                            !        & +a(iix,iiy)*T0x(iix,iiy)**2+b(iix,iiy)*T0y(iix,iiy)**2 &
                            !        & -2*c(iix,iiy)*T0x(iix,iiy)*T0y(iix,iiy) )

                            !sigx=min(1.0,abs(a(iix,iiy)*((px1+px2)/2+T0x(iix,iiy))-c(iix,iiy)*(T0y(iix,iiy)+(py1+py2)/2))/Htau*2)
                            !sigy=min(1.0,abs(b(iix,iiy)*((py1+py2)/2+T0y(iix,iiy))-c(iix,iiy)*(T0x(iix,iiy)+(px1+px2)/2))/Htau*2)
                            
                            !sigx=1!abs(a(iix,iiy)*((px1+px2)/2+T0x(iix,iiy))-c(iix,iiy)*(T0y(iix,iiy)+(py1+py2)/2))/Htau*2
                            !sigy=1!abs(b(iix,iiy)*((py1+py2)/2+T0y(iix,iiy))-c(iix,iiy)*(T0x(iix,iiy)+(px1+px2)/2))/Htau*2
                            Htau = sqrt( a(iix,iiy)*(T0x(iix,iiy)*tau(iix,iiy)+T0(iix,iiy)*(px1+px2)/2)**2 &
                                     & + b(iix,iiy)*(T0y(iix,iiy)*tau(iix,iiy)+T0(iix,iiy)*(py1+py2)/2)**2 &
                                     & - 2*c(iix,iiy)*(T0x(iix,iiy)*tau(iix,iiy)+T0(iix,iiy)*(px1+px2)/2)  &
                                     &               *(T0y(iix,iiy)*tau(iix,iiy)+T0(iix,iiy)*(py1+py2)/2) )

                            tpT=coe*(fun(iix,iiy)-Htau) + coe*(sigx*(px2-px1)/2+sigy*(py2-py1)/2)+tau(iix,iiy);
                                                        
                            tau(iix,iiy)=tpT

                            

                        end if

                    end do
                end do
                
                ! deal with boundary
                do iiy=2,ny-1
                    tau(1,iiy) = max(2*tau(2,iiy)-tau(3,iiy),tau(3,iiy))
                    tau(nx,iiy) = max(2*tau(nx-1,iiy)-tau(nx-2,iiy),tau(nx-2,iiy))
                end do
                do iix=2,nx-1
                    tau(iix,1) = max(2*tau(iix,2)-tau(iix,3),tau(iix,3))
                    tau(iix,ny) = max(2*tau(iix,ny-1)-tau(iix,ny-2),tau(iix,ny-2))
                end do
                
            end do
        end do

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(tau(iix,iiy)-tau_old(iix,iiy))*dx*dy
                Linf_dif=max(Linf_dif,abs(tau(iix,iiy)-tau_old(iix,iiy)))               
            end do
        end do

        do iix=3,nx-2
            do iiy=3,ny-2
                L1_err=L1_err+abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy))
                Linf_err=max(Linf_err,abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy)))            
            end do
        end do

        L1_err = L1_err/(nx-4)/(ny-4)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
        !    write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if
        !write(*,'(a,f10.7,a,f10.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        !write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau*T0

end subroutine

! eikonal equation on a sphere
subroutine FSM_UW_PS_lonlat_2d(xx_deg,yy_deg,nx,ny,spha,sphb,sphc,T,fun,x0_deg,y0_deg,u)
    integer :: nx,ny
    double precision :: xx_deg(nx), yy_deg(ny), x0_deg, y0_deg
    double precision :: xx(nx),yy(ny),spha(nx,ny),sphb(nx,ny),sphc(nx,ny),T(nx,ny),fun(nx,ny),x0,y0,ischange(nx,ny),u(nx,ny)

    double precision :: a(nx,ny),b(nx,ny),c(nx,ny)
    double precision :: dx,dy,T0(nx,ny),T0x(nx,ny),T0y(nx,ny),a0,b0,c0,fun0
    double precision :: tau(nx,ny)
    integer :: iix,iiy,iter,i_sweep, x_id1,x_id2,y_id1,y_id2,xdirec,ydirec
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nx,ny)
    integer,parameter :: MaxIter = 2000
    double precision,parameter :: tol=(10.0)**(-4),eps=(10.0)**(-14)
    double precision,parameter :: R_earth = 6371.0

    integer :: i_case,i_cand,i_solution,count_cand
    double precision :: px,qx,py,qy
    double precision :: eqn_a,eqn_b,eqn_c,Delta,tmp_tau,dis,tmp_T0
    double precision :: characteristic_x, characteristic_y, Tx, Ty
    logical :: is_causility
    double precision :: canditate(50)

    double precision,parameter :: pi=3.14159265358979323846264338327950288

!f2py   intent(in) :: xx_deg, yy_deg
!f2py   intent(hide), depend(xx_deg) :: nx = shape(xx, 0)
!f2py   intent(hide), depend(yy_deg) :: ny = shape(yy, 0)
!f2py   intent(in) :: spha, sphb, sphc, u, fun
!f2py   intent(in) :: x0_deg, y0_deg
!f2py   intent(out) :: T

    ! ------------------------ convert degree to radian ----------------------
    ! do iix=1,nx
    !     xx(iix) = xx_deg(iix)/180.0*pi
    ! end do
    ! do iiy=1,ny
    !     yy(iiy) = yy_deg(iiy)/180.0*pi
    ! end do
    xx = xx_deg/180.0*pi
    yy = yy_deg/180.0*pi
    x0 = x0_deg/180.0*pi
    y0 = y0_deg/180.0*pi

    ! ------------------------ create mesh ------------------------ 
    dx=xx(2)-xx(1); dy=yy(2)-yy(1);

    ! ------------------------ convert a, b, c ----------------------
    ! do iix=1,nx
    !     do iiy=1,ny
    !         a(iix,iiy) = spha(iix,iiy)/(R_earth**2 * cos(yy(iiy))**2)
    !         b(iix,iiy) = sphb(iix,iiy)/(R_earth**2)
    !         c(iix,iiy) = sphc(iix,iiy)/(R_earth**2 * cos(yy(iiy)))
    !     end do
    ! end do
    do iiy = 1,ny
        a(:,iiy) = spha(:,iiy)/(R_earth**2 * cos(yy(iiy))**2)
        b(:,iiy) = sphb(:,iiy)/(R_earth**2)
        c(:,iiy) = sphc(:,iiy)/(R_earth**2 * cos(yy(iiy)))
    end do

    

    ! ------------------------ create T0 ------------------------ 

    ! Parameter discretization at the seismic source.
    idx0=floor((x0-xx(1))/dx+1); idy0=floor((y0-yy(1))/dy+1); 
    r1 = min(1.0,(x0-xx(idx0))/dx); r2 = min(1.0,(y0-yy(idy0))/dy); 

    a0=(1-r1)*(1-r2)*a(idx0,idy0)+(1-r1)*r2*a(idx0,idy0+1) &
    & +r1*(1-r2)*a(idx0+1,idy0)+r1*r2*a(idx0+1,idy0+1) 

    b0=(1-r1)*(1-r2)*b(idx0,idy0)+(1-r1)*r2*b(idx0,idy0+1) &
    & +r1*(1-r2)*b(idx0+1,idy0)+r1*r2*b(idx0+1,idy0+1) 

    c0=(1-r1)*(1-r2)*c(idx0,idy0)+(1-r1)*r2*c(idx0,idy0+1) &
    & +r1*(1-r2)*c(idx0+1,idy0)+r1*r2*c(idx0+1,idy0+1) 

    fun0=(1-r1)*(1-r2)*fun(idx0,idy0)+(1-r1)*r2*fun(idx0,idy0+1) &
    & +r1*(1-r2)*fun(idx0+1,idy0)+r1*r2*fun(idx0+1,idy0+1) 

    ! solve T0 = sqrt((b0(x-x0)^2+a0(y-y0)^2+2c0(x-x0)(y-y0))/(a0b0-c0^2))
    ! write(*,*) fun0,b0/(a0*b0-c0**2),a0/(a0*b0-c0**2),2*c0/(a0*b0-c0**2)
    do iix=1,nx
        do iiy=1,ny
            tmp_T0 = sin(yy(iiy))*sin(y0)+ cos(yy(iiy))*cos(y0)*cos(xx(iix)-x0)
            T0(iix,iiy) = acos(tmp_T0) * R_earth * fun0

            if (T0(iix,iiy) .eq. 0) then
                T0x(iix,iiy) = 0
                T0y(iix,iiy) = 0
            else
                T0x(iix,iiy) = - 1.0 / sqrt(1-tmp_T0**2) * cos(yy(iiy)) *cos(y0) * (-1.0) * sin(xx(iix)-x0) * R_earth * fun0
                T0y(iix,iiy) = - 1.0 / sqrt(1-tmp_T0**2) &
                             & * (cos(yy(iiy))*sin(y0) - sin(yy(iiy))*cos(y0)*cos(xx(iix)-x0)) * R_earth * fun0
            end if
            
            if ( abs((xx(iix)-x0)/dx)<=1 .and. abs((yy(iiy)-y0)/dy)<=1) then
                !Several points around the seismic source are directly assumed to have a constant velocity structure, for which an analytical solution is given, i.e., equal to T0.
                tau(iix,iiy) = 1  
                ischange(iix,iiy)=0
                if (iix==1 .or. iix==nx .or. iiy==1 .or. iiy==ny) then
                    write(*,*) 'source on the boundary, mesh error'
                    write(*,*) x0, y0
                    stop
                end if

            else
                tau(iix,iiy) = 10
                ischange(iix,iiy)=1
            end if


        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do i_sweep = 1,4
            select case(i_sweep)
                case (1)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (2)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1
                case (3)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (4)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1    
            end select


            ! iter 1 x: 1 -> Nx, y: 1 -> Ny
            ! iter 2 x: 1 -> Nx, y: Ny -> 1 
            ! iter 3 x: Nx -> 1, y: 1 -> Ny 
            ! iter 4 x: Nx -> 1, y: Ny -> 1
            !DIR$ SIMD
            do iix=x_id1,x_id2,xdirec
                do iiy=y_id1,y_id2,ydirec
                    if(ischange(iix,iiy)==1) then
    !--------------------------------------- calculate stencil start
                        count_cand = 1
                        do i_case=1,4       
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)    !  -1, -1
                                if (iix == 1 .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)     
                            case (2)    !  -1, +1
                                if (iix == 1 .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1)   
                            case (3)    !  +1, -1
                                if (iix == Nx .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx = +T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                            case (4)    !  +1, +1
                                if (iix == Nx .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                            end select
                            

                            eqn_a = a(iix,iiy)*px**2 + b(iix,iiy)*py**2 -2*c(iix,iiy)*px*py
                            eqn_b = 2*a(iix,iiy)*px*qx + 2*b(iix,iiy)*py*qy - 2*c(iix,iiy)*(px*qy+py*qx)
                            eqn_c = a(iix,iiy)*qx**2 + b(iix,iiy)*qy**2 -2*c(iix,iiy)*qx*qy-fun(iix,iiy)**2
                            Delta = eqn_b**2-4*eqn_a*eqn_c
                            

                            if (Delta>=0) then  !   one or two solutions
                                ! check the causality condition: the characteristic passing through (iix,iiy) is in between used two sides
                                do i_solution=1,2      
                                    select case (i_solution)
                                    case(1)
                                        tmp_tau = (-eqn_b + sqrt(Delta))/(2*eqn_a)
                                    case(2)
                                        tmp_tau = (-eqn_b - sqrt(Delta))/(2*eqn_a)
                                    end select

                                    ! characteristic is
                                    Tx = px*tmp_tau+qx
                                    Ty = py*tmp_tau+qy
                                    characteristic_x = a(iix,iiy)*Tx - c(iix,iiy)*Ty 
                                    characteristic_y = b(iix,iiy)*Ty - c(iix,iiy)*Tx 

                                    is_causility = .false.
                                    ! check the causality condition: 
                                    select case (i_case)
                                    case (1)    
                                        if (characteristic_x >= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (2)   
                                        if (characteristic_x >= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    case (3)    
                                        if (characteristic_x <= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (4)   
                                        if (characteristic_x <= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    end select

                                    ! if satisfying the causility condition, retain it as a canditate solution
                                    if(is_causility) then
                                        canditate(count_cand) = tmp_tau
                                        count_cand = count_cand + 1
                                    end if


                                end do
                            end if

                        end do

                        do i_case=1,4  
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; 
                                ! (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)  
                                if (iix == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy) &
                                & .and. tmp_tau > tau(iix-1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy) &
                                & .and. tmp_tau > tau(iix-1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (2)    
                                if (iix == Nx) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy) &
                                & .and. tmp_tau > tau(iix+1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy) &
                                & .and. tmp_tau > tau(iix+1,iiy) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if

                            case (3)  
                                if (iiy == 1) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1) &
                                & .and. tmp_tau > tau(iix,iiy-1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1) &
                                & .and. tmp_tau > tau(iix,iiy-1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (4)   
                                if (iiy == Ny) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1) &
                                & .and. tmp_tau > tau(iix,iiy+1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1) &
                                & .and. tmp_tau > tau(iix,iiy+1) * 0.5) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            end select

                        end do

                        do i_cand=1,count_cand-1
                            tau(iix,iiy) = min(tau(iix,iiy),canditate(i_cand))
                        end do
    !--------------------------------------- calculate stencil end


                    end if
                end do
            end do
                
                
                
        end do

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy)
                Linf_dif=max(Linf_dif,abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy))               
            end do
        end do

        do iix=1,nx
            do iiy=1,ny
                L1_err=L1_err+abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy))
                Linf_err=max(Linf_err,abs(tau(iix,iiy)*T0(iix,iiy)-u(iix,iiy)))            
            end do
        end do
        L1_dif = L1_dif/(nx)/(ny)
        L1_err = L1_err/(nx)/(ny)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            ! write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
            ! write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if
        ! write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        ! write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau*T0

end subroutine

! Ajoint field, T: traveltime field, Ta: Pn field
! xrec, yrec: locations of receivers
! sourceAdj: traveltime difference, synthetic - observation
! nr: numbers of receivers
subroutine FSM_O1_JSE_2d(xx,yy,nx,ny,a,b,c,T,Ta,xrec,yrec,sourceAdj,nr)
    ! solving \nabla \cdot ( P_n (-\nabla T) M ) = \sum fun*delta(x-x0)   ! only one source
    ! M = [a -c; -c b]
    ! In another word, solveing  (a_new P_n)_x + (b_new P_n)_y = \sum fun*delta(x-x0)
    ! Here a_new= -a T_x + c T_y,  b_new = - b T_y + c T_x 

    integer :: nx,ny,nr
    double precision :: xx(nx),yy(ny),dx,dy,T(nx,ny),sourceAdj(nr),Ta(nx,ny),xrec(nr),yrec(nr),u(nx,ny)
    double precision :: a(nx,ny),b(nx,ny),c(nx,ny)
    double precision :: tpTa,OldTa(nx,ny),d,e
    double precision :: a1,a1m(nx,ny),a1p(nx,ny),a2,a2m(nx,ny),a2p(nx,ny),delta(nx,ny)
    double precision :: b1,b1m(nx,ny),b1p(nx,ny),b2,b2m(nx,ny),b2p(nx,ny)
    integer :: iix,iiy,iter,xdirec,ydirec,idx0,idy0,ir
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    integer,parameter :: MaxIter = 1000
    double precision,parameter :: tol=(10.0)**(-9),eps=(10.0)**(-14)
!f2py   intent(in) :: xx, yy
!f2py   intent(hide) :: nx = shape(xx, 0)
!f2py   intent(hide) :: ny = shape(yy, 0)
!f2py   intent(in) :: a, b, c, T
!f2py   intent(in) :: xrec, yrec, sourceAdj
!f2py   intent(hide) :: nr = shape(xrec, 0)
!f2py   intent(out) :: Ta
    dx=xx(2)-xx(1); dy=yy(2)-yy(1)

    ! create delta source
    do iix=1,nx
        do iiy=1,ny
            delta(iix,iiy)=0
        end do
    end do

    do ir=1,nr
        idx0=floor((xrec(ir)-xx(1))/dx+1); idy0=floor((yrec(ir)-yy(1))/dy+1); 
        r1 = min(1.0,(xrec(ir)-xx(idx0))/dx); r2 = min(1.0,(yrec(ir)-yy(idy0))/dy); 

        delta(idx0,idy0) = delta(idx0,idy0) + sourceAdj(ir)*(1-r1)*(1-r2)/(dx*dy)
        delta(idx0,idy0+1) = delta(idx0,idy0+1) + sourceAdj(ir)*(1-r1)*r2/(dx*dy)
        delta(idx0+1,idy0) = delta(idx0+1,idy0) + sourceAdj(ir)*r1*(1-r2)/(dx*dy)
        delta(idx0+1,idy0+1) = delta(idx0+1,idy0+1) + sourceAdj(ir)*r1*r2/(dx*dy)
    end do

    ! Calculate parameters in advance.
    do iix=2,nx-1
        do iiy=2,ny-1
            a1 = - (T(iix,iiy)-T(iix-1,iiy))/dx * (a(iix,iiy)+a(iix-1,iiy))/2 &
               & + (T(iix,iiy+1)-T(iix,iiy-1)+T(iix-1,iiy+1)-T(iix-1,iiy-1))/(4*dy) &
               &  *(c(iix,iiy)+c(iix-1,iiy))/2
            a1m(iix,iiy) = (a1-abs(a1))/2; a1p(iix,iiy) = (a1+abs(a1))/2;
            a2 = - (T(iix+1,iiy)-T(iix,iiy))/dx * (a(iix+1,iiy)+a(iix,iiy))/2 &
               & + (T(iix+1,iiy+1)-T(iix+1,iiy-1)+T(iix,iiy+1)-T(iix,iiy-1))/(4*dy) &
               &  *(c(iix+1,iiy)+c(iix,iiy))/2
            a2m(iix,iiy) = (a2-abs(a2))/2; a2p(iix,iiy) = (a2+abs(a2))/2;

            b1 = - (T(iix,iiy)-T(iix,iiy-1))/dy * (b(iix,iiy)+b(iix,iiy-1))/2 &
               & + (T(iix+1,iiy)-T(iix-1,iiy)+T(iix+1,iiy-1)-T(iix-1,iiy-1))/(4*dx) &
               &  *(c(iix,iiy)+c(iix,iiy-1))/2
            b1m(iix,iiy) = (b1-abs(b1))/2; b1p(iix,iiy) = (b1+abs(b1))/2;
            b2 = - (T(iix,iiy+1)-T(iix,iiy))/dy * (b(iix,iiy+1)+b(iix,iiy))/2 &
               & + (T(iix+1,iiy+1)-T(iix-1,iiy+1)+T(iix+1,iiy)-T(iix-1,iiy))/(4*dx) &
               &  *(c(iix,iiy+1)+c(iix,iiy))/2
            b2m(iix,iiy) = (b2-abs(b2))/2; b2p(iix,iiy) = (b2+abs(b2))/2;

            Ta(iix,iiy) = 100
        end do
    end do

    ! Initialize the boundary of Ta to 0 and keep it constant.
    do iix=1,nx
        Ta(iix,1) = 0
        Ta(iix,ny) = 0
    end do
    do iiy=1,ny
        Ta(1,iiy) = 0
        Ta(nx,iiy) = 0
    end do

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        OldTa = Ta
        do xdirec = -1,1,2
            do ydirec = -1,1,2
                ! x: 2 <-> Nx-1, y: 2 <-> Ny-1
                
                do iix=nint(0.5+Nx/2.0+(Nx/2.0-1.5)*xdirec),nint(0.5+Nx/2.0+(-Nx/2.0+1.5)*xdirec),-xdirec
                    do iiy=nint(0.5+Ny/2.0+(Ny/2.0-1.5)*ydirec),nint(0.5+Ny/2.0+(-Ny/2.0+1.5)*ydirec),-ydirec

                        d = (a2p(iix,iiy)-a1m(iix,iiy))/dx + (b2p(iix,iiy)-b1m(iix,iiy))/dy

                        if (abs(d)<eps) then
                            Ta(iix,iiy) = 0
                        else
                            e = (Ta(iix-1,iiy)*a1p(iix,iiy)-Ta(iix+1,iiy)*a2m(iix,iiy))/dx &
                            & + (Ta(iix,iiy-1)*b1p(iix,iiy)-Ta(iix,iiy+1)*b2m(iix,iiy))/dy
                            tpTa = (delta(iix,iiy)+e)/d

                            Ta(iix,iiy) = tpTa
                        end if

                    end do
                end do
            end do
        end do

        ! Assess changes.

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(Ta(iix,iiy)-OldTa(iix,iiy))
                Linf_dif=max(Linf_dif,abs(Ta(iix,iiy)-OldTa(iix,iiy)))               
            end do
        end do
        L1_dif = L1_dif/(nx*ny)

        !do iix=1,nx
        !    do iiy=1,ny
        !        L1_err=L1_err+abs(Ta(iix,iiy)-u(iix,iiy))
        !        Linf_err=max(Linf_err,abs(Ta(iix,iiy)-u(iix,iiy)))               
        !    end do
        !end do
        !L1_err = L1_err/(nx*ny)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            ! write(*,*) 'iter ',iter,', Ta is steady'
            exit
        ! else
        !    write(*,*) 'iter ',iter,', Ta is changing, continue ... '
        end if

        ! if (iter==MaxIter) then    
        !     write(*,*) 'iter ',iter,', max iteration steps'
        ! end if
        !write(*,'(a,f10.7,a,f10.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        !write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err


    end do


end subroutine

! adjoint field on a sphere
subroutine FSM_O1_JSE_lonlat_2d(xx_deg,yy_deg,nx,ny,spha,sphb,sphc,T,Ta,xrec_deg,yrec_deg,sourceAdj,nr)
    ! solving \nabla \cdot ( P_n (-\nabla T) M ) = \sum fun*delta(x-x0)   ! only one source
    ! M = [a -c; -c b]
    ! In another word, solveing  (a_new P_n)_x + (b_new P_n)_y = \sum fun*delta(x-x0)
    ! Here a_new= -a T_x + c T_y,  b_new = - b T_y + c T_x 

    integer :: nx,ny,nr
    double precision :: xx_deg(nx),yy_deg(ny),xrec_deg(nr),yrec_deg(nr)
    double precision :: xx(nx),yy(ny),dx,dy,T(nx,ny),sourceAdj(nr),Ta(nx,ny),xrec(nr),yrec(nr),u(nx,ny)
    double precision :: a(nx,ny),b(nx,ny),c(nx,ny)

    double precision :: spha(nx,ny),sphb(nx,ny),sphc(nx,ny)
    double precision :: tpTa,OldTa(nx,ny),d,e
    double precision :: a1,a1m(nx,ny),a1p(nx,ny),a2,a2m(nx,ny),a2p(nx,ny),delta(nx,ny)
    double precision :: b1,b1m(nx,ny),b1p(nx,ny),b2,b2m(nx,ny),b2p(nx,ny)
    integer :: iix,iiy,iter,xdirec,ydirec,idx0,idy0,ir
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    integer,parameter :: MaxIter = 1000
    double precision,parameter :: tol=(10.0)**(-9),eps=(10.0)**(-14)

    double precision,parameter :: R_earth = 6371.0
    double precision,parameter :: pi=3.14159265358979323846264338327950288

!f2py   intent(in) :: xx_deg, yy_deg
!f2py   intent(hide) :: nx = shape(xx_deg, 0)
!f2py   intent(hide) :: ny = shape(yy_deg, 0)
!f2py   intent(in) :: spha, sphb, sphc, T
!f2py   intent(in) :: xrec_deg, yrec_deg, sourceAdj
!f2py   intent(hide) :: nr = shape(xrec_deg, 0)
!f2py   intent(out) :: Ta

    ! convert degree to radian
    do iix=1,nx
        xx(iix) = xx_deg(iix)/180.0*pi
    end do
    do iiy=1,ny
        yy(iiy) = yy_deg(iiy)/180.0*pi
    end do
    do ir = 1,nr
        xrec(ir) = xrec_deg(ir)/180.0*pi
        yrec(ir) = yrec_deg(ir)/180.0*pi
    end do

    dx=xx(2)-xx(1); dy=yy(2)-yy(1)

    ! create delta source
    do iix=1,nx
        do iiy=1,ny
            delta(iix,iiy)=0
        end do
    end do

    do ir=1,nr
        idx0=floor((xrec(ir)-xx(1))/dx+1); idy0=floor((yrec(ir)-yy(1))/dy+1); 
        r1 = min(1.0,(xrec(ir)-xx(idx0))/dx); r2 = min(1.0,(yrec(ir)-yy(idy0))/dy); 

        delta(idx0,idy0) = delta(idx0,idy0) + sourceAdj(ir)*(1-r1)*(1-r2)/(dx*R_earth*cos(yrec(ir))*dy*R_earth)
        delta(idx0,idy0+1) = delta(idx0,idy0+1) + sourceAdj(ir)*(1-r1)*r2/(dx*R_earth*cos(yrec(ir))*dy*R_earth)
        delta(idx0+1,idy0) = delta(idx0+1,idy0) + sourceAdj(ir)*r1*(1-r2)/(dx*R_earth*cos(yrec(ir))*dy*R_earth)
        delta(idx0+1,idy0+1) = delta(idx0+1,idy0+1) + sourceAdj(ir)*r1*r2/(dx*R_earth*cos(yrec(ir))*dy*R_earth)
    end do

    do iix=1,nx
        do iiy=1,ny
            a(iix,iiy) = spha(iix,iiy)/(R_earth**2 * cos(yy(iiy))**2)
            b(iix,iiy) = sphb(iix,iiy)/(R_earth**2)
            c(iix,iiy) = sphc(iix,iiy)/(R_earth**2 * cos(yy(iiy)))
        end do
    end do


    ! Calculate parameters in advance.
    do iix=2,nx-1
        do iiy=2,ny-1
            a1 = - (T(iix,iiy)-T(iix-1,iiy))/dx * (a(iix,iiy)+a(iix-1,iiy))/2 &
               & + (T(iix,iiy+1)-T(iix,iiy-1)+T(iix-1,iiy+1)-T(iix-1,iiy-1))/(4*dy) &
               &  *(c(iix,iiy)+c(iix-1,iiy))/2
            a1m(iix,iiy) = (a1-abs(a1))/2; a1p(iix,iiy) = (a1+abs(a1))/2;
            a2 = - (T(iix+1,iiy)-T(iix,iiy))/dx * (a(iix+1,iiy)+a(iix,iiy))/2 &
               & + (T(iix+1,iiy+1)-T(iix+1,iiy-1)+T(iix,iiy+1)-T(iix,iiy-1))/(4*dy) &
               &  *(c(iix+1,iiy)+c(iix,iiy))/2
            a2m(iix,iiy) = (a2-abs(a2))/2; a2p(iix,iiy) = (a2+abs(a2))/2;

            b1 = - (T(iix,iiy)-T(iix,iiy-1))/dy * (b(iix,iiy)+b(iix,iiy-1))/2 &
               & + (T(iix+1,iiy)-T(iix-1,iiy)+T(iix+1,iiy-1)-T(iix-1,iiy-1))/(4*dx) &
               &  *(c(iix,iiy)+c(iix,iiy-1))/2
            b1m(iix,iiy) = (b1-abs(b1))/2; b1p(iix,iiy) = (b1+abs(b1))/2;
            b2 = - (T(iix,iiy+1)-T(iix,iiy))/dy * (b(iix,iiy+1)+b(iix,iiy))/2 &
               & + (T(iix+1,iiy+1)-T(iix-1,iiy+1)+T(iix+1,iiy)-T(iix-1,iiy))/(4*dx) &
               &  *(c(iix,iiy+1)+c(iix,iiy))/2
            b2m(iix,iiy) = (b2-abs(b2))/2; b2p(iix,iiy) = (b2+abs(b2))/2;

            Ta(iix,iiy) = 100
        end do
    end do

    ! Initialize the boundary of Ta to 0 and keep it constant.
    do iix=1,nx
        Ta(iix,1) = 0
        Ta(iix,ny) = 0
    end do
    do iiy=1,ny
        Ta(1,iiy) = 0
        Ta(nx,iiy) = 0
    end do

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        OldTa = Ta
        do xdirec = -1,1,2
            do ydirec = -1,1,2
                ! x: 2 <-> Nx-1, y: 2 <-> Ny-1
                
                do iix=nint(0.5+Nx/2.0+(Nx/2.0-1.5)*xdirec),nint(0.5+Nx/2.0+(-Nx/2.0+1.5)*xdirec),-xdirec
                    do iiy=nint(0.5+Ny/2.0+(Ny/2.0-1.5)*ydirec),nint(0.5+Ny/2.0+(-Ny/2.0+1.5)*ydirec),-ydirec

                        d = (a2p(iix,iiy)-a1m(iix,iiy))/dx + (b2p(iix,iiy)-b1m(iix,iiy))/dy

                        if (abs(d)<eps) then
                            Ta(iix,iiy) = 0
                        else
                            e = (Ta(iix-1,iiy)*a1p(iix,iiy)-Ta(iix+1,iiy)*a2m(iix,iiy))/dx &
                            & + (Ta(iix,iiy-1)*b1p(iix,iiy)-Ta(iix,iiy+1)*b2m(iix,iiy))/dy
                            tpTa = (delta(iix,iiy)+e)/d

                            Ta(iix,iiy) = tpTa
                        end if

                    end do
                end do

            end do
        end do

        ! Assess changes.

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(Ta(iix,iiy)-OldTa(iix,iiy))
                Linf_dif=max(Linf_dif,abs(Ta(iix,iiy)-OldTa(iix,iiy)))               
            end do
        end do
        L1_dif = L1_dif/(nx*ny)

        !do iix=1,nx
        !    do iiy=1,ny
        !        L1_err=L1_err+abs(Ta(iix,iiy)-u(iix,iiy))
        !        Linf_err=max(Linf_err,abs(Ta(iix,iiy)-u(iix,iiy)))               
        !    end do
        !end do
        !L1_err = L1_err/(nx*ny)

        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            ! write(*,*) 'iter ',iter,', Ta is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', Ta is changing, continue ... '
        end if

        if (iter==MaxIter) then    
            ! write(*,*) 'iter ',iter,', max iteration steps'
        end if
        !write(*,'(a,f10.7,a,f10.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        !write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err


    end do


end subroutine


subroutine mask(xx,yy,nx,ny,TableAdj,x0,y0)
    integer :: nx,ny
    double precision :: xx(nx),yy(ny),TableAdj(nx,ny),x0,y0,dx,dy
    integer :: iix,iiy
!f2py   intent(in) :: xx, yy
!f2py   intent(hide), depend(xx) :: nx = shape(xx, 0)
!f2py   intent(hide), depend(yy) :: ny = shape(yy, 0)
!f2py   intent(in) :: x0, y0
!f2py   intent(inout) :: TableAdj
    dx=xx(2)-xx(1); dy=yy(2)-yy(1)

    do iix=1,nx
        do iiy=1,ny
            if ( abs(xx(iix)-x0)<1.5*dx .and. abs(yy(iiy)-y0)<1.5*dy ) then
                TableAdj(iix,iiy) = 0
            end if
        end do
    end do
end subroutine


subroutine locate_within_cube(xx, yy, zz, nx, ny, nz, xvel, yvel, zvel, &
    i1, j1, k1, wx, wy, wz, wt, iloc)
    integer :: nx, ny, nz, i1, j1, k1, iloc
    double precision :: xx, yy, zz, wx, wy, wz, wt(8)
    double precision :: xvel(nx), yvel(ny), zvel(nz)

    integer :: i2, j2, k2
    logical :: notfindx, notfindy, notfindz

    iloc = 0
    notfindx = .true.
    do while (notfindx .and. i1 .gt. 0 .and. i1 .lt. nx)
    i2 = i1 + 1
    if (xx .lt. xvel(i1)) then
    i1 = i1 - 1
    else
    if (xx .le. xvel(i2)) then
    notfindx = .false.
    wx = (xx - xvel(i1)) / (xvel(i2) - xvel(i1))
    else
    i1 = i1 + 1
    end if
    end if
    end do

    notfindy = .true.
    do while (notfindy .and. j1 .gt. 0 .and. j1 .lt. ny)
    j2 = j1 + 1
    if (yy .lt. yvel(j1)) then
    j1 = j1 - 1
    else
    if (yy .le. yvel(j2)) then
    notfindy = .false.
    wy = (yy - yvel(j1)) / (yvel(j2) - yvel(j1))
    else
    j1 = j1 + 1
    end if
    end if
    end do

    notfindz = .true.
    do while (notfindz .and. k1 .gt. 0 .and. k1 .lt. nz)
    k2 = k1 + 1
    if (zz .lt. zvel(k1)) then
    k1 = k1 - 1
    else
    if (zz .le. zvel(k2)) then
    notfindz = .false.
    wz = (zz - zvel(k1)) / (zvel(k2) - zvel(k1))
    else
    k1 = k1 + 1
    end if
    end if
    end do

    if ((.not.notfindx) .and. (.not.notfindy) .and. (.not.notfindz)) then
    iloc  = 1
    wt(1) = (1.0 - wx) * (1.0 - wy) * (1.0 - wz)
    wt(2) = (1.0 - wx) * wy * (1.0 - wz)
    wt(3) = wx * wy * (1.0 - wz)
    wt(4) = wx * (1.0 - wy) * (1.0 - wz)
    wt(5) = (1.0 - wx) * (1.0 - wy) * wz
    wt(6) = (1.0 - wx)* wy * wz
    wt(7) = wx * wy * wz
    wt(8) = wx * (1.0 - wy) * wz
    end if

    return
end subroutine locate_within_cube


subroutine inv_grid_iso(xinv, yinv, zinv, nxinv, nyinv, nzinv, nset, adj, nx, ny, nz_3d, gk, vel_x, vel_y, vel_z_3d)
    ! this subroutine cannot deal with multiple scale grid
    integer iset, nset, inside, n, m
    integer nxinv, nyinv, nzinv
    integer i1, j1, k1, i, j, k, nx1, ny1, nz1
    integer nx, ny, nz_3d
    integer istart
    double precision :: wx, wy, wz, wt(8)
    double precision :: xinv(nxinv, nset), yinv(nyinv, nset), zinv(nzinv, nset)
    double precision :: xinvtmp(nxinv), yinvtmp(nyinv), zinvtmp(nzinv)
    double precision :: adj(nx, ny, nz_3d)
    double precision :: gk(nset*nxinv*nyinv*nzinv)
    double precision :: vel_x(nx), vel_y(ny), vel_z_3d(nz_3d), dx, dy, dz
!f2py   intent(in) :: xinv, yinv, zinv, adj, vel_x, vel_y, vel_z_3d
!f2py   intent(out) :: gk
!f2py   intent(hide), depend(xinv)   :: nset=shape(xinv, 1), nxinv=shape(xinv,0)
!f2py   intent(hide), depend(yinv)   :: nyinv=shape(yinv,0)
!f2py   intent(hide), depend(zinv)   :: nzinv=shape(zinv,0)
!f2py   intent(hide), depend(adj) :: nx=shape(adj, 0), ny=shape(adj, 1), nz_3d=shape(adj, 2)
    dx = vel_x(2) - vel_x(1)
    dy = vel_y(2) - vel_y(1)
    dz = vel_z_3d(2) - vel_z_3d(1)
    nx1 = nx - 1
    ny1 = ny - 1
    nz1 = nz_3d - 1
    istart = 0
    gk = 0
    do iset = 1, nset
        xinvtmp = xinv(:, iset)
        yinvtmp = yinv(:, iset)
        zinvtmp = zinv(:, iset)
        i1 = 1
        j1 = 1
        k1 = 1

        do j = 2, ny1
            do i = 2, nx1
                do k = 2, nz1
                    call locate_within_cube(vel_x(i),vel_y(j),vel_z_3d(k),nxinv, nyinv, nzinv,&
                                            xinvtmp,yinvtmp,zinvtmp,i1,j1,k1,wx,wy,wz,wt,inside)

                    if (inside .gt. 0) then
                        do n = 1, 8
                            if (n .eq. 1) then
                                m = nxinv*nyinv*(k1-1)+nxinv*(j1-1)+i1
                            elseif (n .eq. 2) then
                                m = nxinv*nyinv*(k1-1)+nxinv*j1+i1
                            elseif (n .eq. 3) then
                                m = nxinv*nyinv*(k1-1)+nxinv*j1+i1+1
                            elseif (n .eq. 4) then
                                m = nxinv*nyinv*(k1-1)+nxinv*(j1-1)+i1+1
                            elseif (n .eq. 5) then
                                m = nxinv*nyinv*k1+nxinv*(j1-1)+i1
                            elseif (n .eq. 6) then
                                m = nxinv*nyinv*k1+nxinv*j1+i1
                            elseif (n .eq. 7) then
                                m = nxinv*nyinv*k1+nxinv*j1+i1+1
                            elseif (n .eq. 8) then
                                m = nxinv*nyinv*k1+nxinv*(j1-1)+i1+1
                            else
                            end if
                            m = istart + m
                            gk(m) = gk(m) + adj(i,j,k)*dx*dy*dz*wt(n)
                        end do
                    else
                        print *, "I cannot locate the receiver in the domain.", vel_x(i),vel_y(j),vel_z_3d(k)
                        stop
                    end if        
                end do        
            end do
        end do
        istart = istart + nxinv * nyinv * nzinv
    end do
end subroutine inv_grid_iso


subroutine inv2fwd_iso(xinv, yinv, zinv, nxinv, nyinv, nzinv, nset, nx, ny, nz_3d, gk, vel_x, vel_y, vel_z_3d, update)
    ! this subroutine cannot deal with multiple scale grid
    integer iset, nset, inside, n, m
    integer nxinv, nyinv, nzinv
    integer i1, j1, k1, i, j, k
    integer nx, ny, nz_3d
    integer istart
    double precision :: wx, wy, wz, wt(8)
    double precision :: xinv(nxinv, nset), yinv(nyinv, nset), zinv(nzinv, nset)
    double precision :: xinvtmp(nxinv), yinvtmp(nyinv), zinvtmp(nzinv)
    double precision :: gk(nset*nxinv*nyinv*nzinv)
    double precision :: vel_x(nx), vel_y(ny), vel_z_3d(nz_3d)
    double precision :: update(nx, ny, nz_3d)
    double precision :: val
!f2py   intent(in) :: xinv, yinv, zinv, vel_x, vel_y, vel_z_3d, gk
!f2py   intent(out) :: update
!f2py   intent(hide), depend(xinv)   :: nset=shape(xinv, 1), nxinv=shape(xinv,0)
!f2py   intent(hide), depend(yinv)   :: nyinv=shape(yinv,0)
!f2py   intent(hide), depend(zinv)   :: nzinv=shape(zinv,0)
!f2py   intent(hide), depend(vel_x) :: nx=shape(vel_x, 0)
!f2py   intent(hide), depend(vel_y) :: ny=shape(vel_y, 0)
!f2py   intent(hide), depend(vel_z_3d) :: nz_3d=shape(vel_z_3d, 0)

    istart = 0
    update = 0

    do iset = 1, nset
        xinvtmp = xinv(:, iset)
        yinvtmp = yinv(:, iset)
        zinvtmp = zinv(:, iset)
        i1 = 1
        j1 = 1
        k1 = 1

        do j = 1, ny
            do i = 1, nx
                do k = 1, nz_3d
                    call locate_within_cube(vel_x(i),vel_y(j),vel_z_3d(k),nxinv, nyinv, nzinv,&
                                            xinvtmp,yinvtmp,zinvtmp,i1,j1,k1,wx,wy,wz,wt,inside)

                    if (inside .gt. 0) then
                        val = 0.0
                        do n = 1, 8
                            if (n .eq. 1) then
                                m = nxinv*nyinv*(k1-1)+nxinv*(j1-1)+i1
                            elseif (n .eq. 2) then
                                m = nxinv*nyinv*(k1-1)+nxinv*j1+i1
                            elseif (n .eq. 3) then
                                m = nxinv*nyinv*(k1-1)+nxinv*j1+i1+1
                            elseif (n .eq. 4) then
                                m = nxinv*nyinv*(k1-1)+nxinv*(j1-1)+i1+1
                            elseif (n .eq. 5) then
                                m = nxinv*nyinv*k1+nxinv*(j1-1)+i1
                            elseif (n .eq. 6) then
                                m = nxinv*nyinv*k1+nxinv*j1+i1
                            elseif (n .eq. 7) then
                                m = nxinv*nyinv*k1+nxinv*j1+i1+1
                            elseif (n .eq. 8) then
                                m = nxinv*nyinv*k1+nxinv*(j1-1)+i1+1
                            else
                            end if
                            m = istart + m
                            val = val + wt(n) * gk(m)
                        end do
                    else
                        print *, "I cannot locate the receiver in the domain.", vel_x(i),vel_y(j),vel_z_3d(k)
                        stop
                    end if
                    update(i, j, k) = update(i, j, k) + val        
                end do        
            end do
        end do
        istart = istart + nxinv * nyinv * nzinv
    end do
end subroutine inv2fwd_iso

