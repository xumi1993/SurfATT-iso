subroutine depthkernelTI_mpi(vel,nx,ny,nz,igrid,istart,iend,&
                      kmaxRc, tRc, cgRc, depz, dcR_dA, dcR_dL)
  use surfker, only: get_vprho
  implicit none
  integer, intent(in) :: nx, ny, nz, istart, iend
  double precision, dimension(nx,ny,nz), intent(in) :: vel
  integer, intent(in) :: kmaxRc
  double precision, dimension(kmaxRc), intent(in) :: tRc
  double precision, dimension(kmaxRc,nx,ny), intent(in) :: cgRc!svel
  double precision, dimension(nz), intent(in) :: depz
  integer, dimension(nx*ny, 2), intent(in) :: igrid
  double precision, dimension(kmaxRc, nx, ny, nz), intent(out) :: dcR_dL, dcR_dA
  integer :: mmax, iflsph, mode, ii, ix, iy, i
  real, dimension(nz) :: vpz,vsz,rhoz
  integer, parameter :: NL=200, NP=60
  real, dimension(NP) :: cp_in, t_in
  real, dimension(NL) :: TA_in,TC_in,TL_in,TN_in,TF_in,TRho_in,rthk,qp,qs,etap,etas,frefp,frefs
  real, dimension(NP,NL) :: dcdah,dcdn,dcdbv

  mmax=nz
  iflsph=1
  mode=1
  dcR_dA = 0.
  dcR_dL = 0.

  do ii = istart, iend
    ix = igrid(ii,1)
    iy = igrid(ii,2)
    vsz(1:nz)=sngl(vel(ix,iy,1:nz))
    call get_vprho(vsz, nz, vpz, rhoz)
    do i = 1, mmax
      TA_in(i)=rhoz(i)*vpz(i)**2
      TC_in(i)=TA_in(i)
      TL_in(i)=rhoz(i)*vsz(i)**2
      TN_in(i)=TL_in(i)
      TF_in(i)=1.0*(TA_in(i) - 2 * TL_in(i))
      TRho_in(i)=rhoz(i)
    enddo
    do i = 1, mmax-1
      rthk(i) = sngl(depz(i+1) - depz(i))
    enddo
    rthk(mmax) = 0.
    qp(1:mmax)=150.0
    qs(1:mmax)=50.0
    etap(1:mmax)=0.00
    etas(1:mmax)=0.00
    frefp(1:mmax)=1.00
    frefs(1:mmax)=1.00
    cp_in(1:kmaxRc)=sngl(cgRc(1:kmaxRc, ix, iy))
    t_in(1:kmaxRc)=sngl(tRc(1:kmaxRc))

    call tregn96(mmax, rthk, TA_in, TC_in, TF_in,&
                 TL_in, TN_in, TRho_in, &
                 qp, qs, etap, etas, frefp, frefs,  &
                 kmaxRc, t_in, cp_in,&
                 dcdah, dcdbv, dcdn)
    do i=1,kmaxRc  ! period
      dcR_dA(i, ix, iy, :) = 0.5/(rhoz*vpz)*dcdah(i,1:nz) -&
          TF_in(1:nz)/((TA_in(1:nz)-2.0*TL_in(1:nz))**2)*dcdn(i,1:nz)
      dcR_dL(i, ix, iy, :) = 0.5/(rhoz*vsz)*dcdbv(i, 1:nz) +&
          2.0*TF_in(1:nz)/((TA_in(1:nz)-2.0*TL_in(1:nz))**2)*dcdn(i,1:nz)
    enddo
    ! dcR_dA(:,:,:,nz) = 0.
    ! dcR_dL(:,:,:,nz) = 0.
  enddo
    
end subroutine depthkernelTI_mpi
! use tregn96 from cps to calculate the dcdL and dcdA
! subroutine depthkernelTI(nx, ny, nz, vel, iwave, igr, kmaxRc, tRc, depz, minthk, dcR_dA, dcR_dL, pvRc)
subroutine depthkernelTI(nx, ny, nz, vel, iwave, igr, kmaxRc, tRc, depz, minthk, dcR_dA, dcR_dL)
    use constants
    implicit none

    integer nx, ny, nz
    real vel(nx, ny, nz)
    integer iwave,igr
    real minthk
    real depz(nz)
    integer kmaxRc
    real(kind=dp) tRc(kmaxRc)
    ! output
    real(kind=dp) pvRc(kmaxRc, nx, ny)
    ! parameter list
    real vpz(nz),vsz(nz),rhoz(nz)
    integer mmax,iflsph,mode,rmax
    integer ii,jj,k,i,j,jjj
    integer, parameter:: NL=200
    integer, parameter:: NP=60

    real(kind=dp) cgRc(NP)
    real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    ! for tregn96
    real t_in(kmaxRc), cp_in(kmaxRc)
    real TA_in(NL), TC_in(NL), TF_in(NL)
    real TL_in(NL), TN_in(NL), TRho_in(NL)
    real qp(NL), qs(NL), etap(NL)
    real etas(NL), frefp(NL), frefs(NL)

    real(kind=cr) dcdah(NP,NL),dcdn(NP,NL)
    real(kind=cr) dcdbv(NP,NL)
    real(kind=dp) dcR_dL(kmaxRc, nx, ny, nz), dcR_dA(kmaxRc, nx, ny, nz)
    integer nsublay(NL), post ! unused

!f2py intent(in) :: vel, iwave, igr, tRc, depz, minthk
!f2py intent(out) :: dcR_dA
!f2py intent(out) :: dcR_dL
!f2py intent(out) :: pvRc
!f2py intent(hide), depend(vel) :: nx=shape(vel, 0), ny=shape(vel, 1), nz=shape(vel, 2)
!f2py intent(hide), depend(tRc) :: kmaxRc=shape(tRc, 0)

    mmax=nz
    iflsph=1
    mode=1
    ! pvRc=0.0

    do jj=1,ny
        do ii=1,nx
            ! post=ii+(jj-1)*nx ! unused
            vsz(1:nz)=vel(ii,jj,1:nz)
            ! some other emperical relationship maybe better,
            do k=1,nz
                vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
                0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
                rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
                0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
                0.000106*vpz(k)**5
            enddo
            ! change from refineGrid2LayerMdl into refineLayerMdl
            ! call refineGrid2LayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            ! rvp, rvs, rrho, rthk)
            call refineLayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            rvp, rvs, rrho, rthk, nsublay)

            call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, kmaxRc, &
            tRc, cgRc)
            ! pvRc(1:kmaxRc, ii, jj)=cgRc(1:kmaxRc)
            !print*,cgRc(1:kmaxRc)
            !------------------------------------------------------------------!
            do i = 1, rmax
                TA_in(i)=rrho(i)*rvp(i)**2
                TC_in(i)=TA_in(i)
                TL_in(i)=rrho(i)*rvs(i)**2
                TN_in(i)=TL_in(i)
                TF_in(i)=1.0*(TA_in(i) - 2 * TL_in(i))
                TRho_in(i)=rrho(i)
            enddo
            qp(1:rmax)=150.0
            qs(1:rmax)=50.0
            etap(1:rmax)=0.00
            etas(1:rmax)=0.00
            frefp(1:rmax)=1.00
            frefs(1:rmax)=1.00

            cp_in(1:kmaxRc)=sngl(cgRc(1:kmaxRc))
            t_in(1:kmaxRc)=sngl(tRc(1:kmaxRc))

            ! ! write(6, *)'tregn96'
            call tregn96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
            qp, qs, etap, etas, frefp, frefs,  &
            kmaxRc, t_in, cp_in(1:kmaxRc),&
            dcdah, dcdbv, dcdn)
            !
            ! ! write(*,*)"nsublay:", nsublay(1:nz)
            do i=1,kmaxRc  ! period
                k=0
                do j=1,nz-1                ! inversion layer
                    do jjj=1,nsublay(j)    ! refined layer k-th in jth inversion layer
                        k=k+1
                        !TODO change refineLayerMdl, use layer thickness same as input
                        dcR_dA(i, ii, jj, j) = 0.5/(rrho(k)*rvp(k))*dcdah(i, k) -&
                            TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                        dcR_dL(i, ii, jj, j) = 0.5/(rrho(k)*rvs(k))*dcdbv(i, k) +&
                            2.0*TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                    enddo
                enddo
            enddo
        enddo
    enddo

end subroutine

subroutine depthkernelTI1d(nz, vel, iwave, igr, kmaxRc, tRc, depz, minthk, dcR_dA, dcR_dL, pvRc)
    use constants
    implicit none

    integer nz
    real vel(nz)
    integer iwave,igr
    real minthk
    real depz(nz)
    integer kmaxRc
    real(kind=dp) tRc(kmaxRc)
    ! output
    real(kind=dp) pvRc(kmaxRc)
    ! parameter list
    real vpz(nz),vsz(nz),rhoz(nz)
    integer mmax,iflsph,mode,rmax
    integer ii,jj,k,i,j,jjj
    integer, parameter:: NL=200
    integer, parameter:: NP=60

    real(kind=dp) cgRc(NP)
    real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    ! for tregn96
    real t_in(kmaxRc), cp_in(kmaxRc)
    real TA_in(NL), TC_in(NL), TF_in(NL)
    real TL_in(NL), TN_in(NL), TRho_in(NL)
    real qp(NL), qs(NL), etap(NL)
    real etas(NL), frefp(NL), frefs(NL)

    real(kind=cr) dcdah(NP,NL),dcdn(NP,NL)
    real(kind=cr) dcdbv(NP,NL)
    real(kind=cr) dcR_dL(kmaxRc, nz), dcR_dA(kmaxRc, nz)
    integer nsublay(NL), post

!f2py intent(in) :: vel, iwave, igr, tRc, depz, minthk
!f2py intent(out) :: dcR_dA, dcR_dL, pvRc
!f2py intent(hide), depend(vel) :: nz=shape(vel, 0)
!f2py intent(hide), depend(tRc) :: kmaxRc=shape(tRc, 0)

    mmax=nz
    iflsph=1
    mode=1
    pvRc=0.0


    vsz(1:nz)=vel(1:nz)
    ! some other emperical relationship maybe better,
    do k=1,nz
        vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
        0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
        rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
        0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
        0.000106*vpz(k)**5
    enddo
    ! change from refineGrid2LayerMdl into refineLayerMdl
    ! call refineGrid2LayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
    ! rvp, rvs, rrho, rthk)
    call refineLayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
    rvp, rvs, rrho, rthk, nsublay)

    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, kmaxRc, &
    tRc, cgRc)
    pvRc(1:kmaxRc)=cgRc(1:kmaxRc)
    !print*,cgRc(1:kmaxRc)
    !------------------------------------------------------------------!
    do i = 1, rmax
        TA_in(i)=rrho(i)*rvp(i)**2
        TC_in(i)=TA_in(i)
        TL_in(i)=rrho(i)*rvs(i)**2
        TN_in(i)=TL_in(i)
        TF_in(i)=1.0*(TA_in(i) - 2 * TL_in(i))
        TRho_in(i)=rrho(i)
    enddo
    qp(1:rmax)=150.0
    qs(1:rmax)=50.0
    etap(1:rmax)=0.00
    etas(1:rmax)=0.00
    frefp(1:rmax)=1.00
    frefs(1:rmax)=1.00

    cp_in(1:kmaxRc)=sngl(cgRc(1:kmaxRc))
    t_in(1:kmaxRc)=sngl(tRc(1:kmaxRc))

    ! ! write(6, *)'tregn96'
    call tregn96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
    qp, qs, etap, etas, frefp, frefs,  &
    kmaxRc, t_in, cp_in(1:kmaxRc),&
    dcdah, dcdbv, dcdn)
    !
    ! ! write(*,*)"nsublay:", nsublay(1:nz)
    do i=1,kmaxRc  ! period
        k=0
        do j=1,nz-1                ! inversion layer
            do jjj=1,nsublay(j)    ! refined layer k-th in jth inversion layer
                k=k+1
                dcR_dA(i,j) = 0.5/(rrho(k)*rvp(k))*dcdah(i, k) - TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                dcR_dL(i,j) = 0.5/(rrho(k)*rvs(k))*dcdbv(i, k) + 2.0*TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
            enddo
        enddo
    enddo
end subroutine


subroutine refineLayerMdl(minthk0,mmax,dep,vp,vs,rho,&
    rmax,rdep,rvp,rvs,rrho,rthk,nsublay)
!!--------------------------------------------------------------------c
!!refine grid based model to layerd based model
!!input:   minthk: is the minimum thickness of the refined layered model
!!         mmax: number of depth grid points in the model
!!         dep, vp, vs, rho: the depth-grid model parameters
!!         rmax: number of layers in the fined layered model
!!         rdep, rvp, rvs, rrho, rthk: the refined layered velocity model
!!
implicit none
integer NL
parameter (NL=200)
integer mmax,rmax
real minthk0
real minthk
real dep(*),vp(*),vs(*),rho(*)
real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
integer nsublay(NL)
real thk,newthk,initdep
integer i,j,k,ngrid

k = 0
initdep = 0.0
do i = 1, mmax-1
thk = dep(i+1)-dep(i)
minthk = thk/minthk0
nsublay(i) = int((thk+1.0e-4)/minthk) + 1
ngrid = nsublay(i)+1
newthk = thk/nsublay(i)
do j = 1, nsublay(i)
k = k + 1
rthk(k) = newthk
rdep(k) = initdep + rthk(k)
initdep = rdep(k)
rvp(k) = vp(i)+(2*j-1)*(vp(i+1)-vp(i))/(2*nsublay(i))
rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*nsublay(i))
rrho(k) = rho(i)+(2*j-1)*(rho(i+1)-rho(i))/(2*nsublay(i))
enddo
enddo
!! half space model
!        k = k + 1
!        rthk(k) = 0.0
!        rvp(k) = vp(mmax)
!        rvs(k) = vs(mmax)
!        rrho(k) = rho(mmax)
!	rdep(k) = dep(mmax)
rmax = k
return
end subroutine

! subroutine fwdsurf1d(vel,nz,iwave,igr,kmaxRc,tRc,depz,minthk,svel)
!     implicit none
  
!     integer nz
!     real vel(nz)
!     real(kind=dp) svel(kmaxRc)
  
!     integer iwave,igr
!     real minthk
!     real depz(nz)
!     integer kmaxRc
!     real(kind=dp) tRc(kmaxRc)
  
!     real vpz(nz),vsz(nz),rhoz(nz)
!     integer mmax,iflsph,mode,rmax
!     integer k
!     integer,parameter::NL=200
!     integer,parameter::NP=100
!     real(kind=dp) cgRc(NP)
!     real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
  
!   !f2py intent(in) :: vel, iwave, igr, tRc, depz, minthk
!   !f2py intent(out) :: svel
!   !f2py intent(hide), depend(vel) :: nz=shape(vel, 0)
!   !f2py intent(hide), depend(tRc) :: kmaxRc=shape(tRc, 0)
    
!     mmax=nz
!     iflsph=1
!     mode=1
  
!     vsz(1:nz)=vel(1:nz)
!     ! some other emperical relationship maybe better, 
!     ! This is from Tomas M.Brocher 2005 BSSA
!     do k=1,nz
!       vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
!         0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
!       rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
!         0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + & 
!         0.000106*vpz(k)**5
!     enddo
  
!     call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax,rdep,&
!       rvp,rvs,rrho,rthk)
!     call surfdisp96(rthk,rvp,rvs,rrho,rmax,iflsph,iwave,mode,igr,kmaxRc,&
!       tRc,cgRc)
!     svel(1:kmaxRc)=cgRc(1:kmaxRc)
  
!   end subroutine fwdsurf1d
  
  ! subroutine refineGrid2LayerMdl(minthk0,mmax,dep,vp,vs,rho,&
  !     rmax,rdep,rvp,rvs,rrho,rthk)
  !   !!--------------------------------------------------------------------c
  !   !!refine grid based model to layerd based model
  !   !!INPUT:   minthk: is the minimum thickness of the refined layered model
  !   !!         mmax: number of depth grid points in the model
  !   !!         dep, vp, vs, rho: the depth-grid model parameters
  !   !!         rmax: number of layers in the fined layered model
  !   !!         rdep, rvp, rvs, rrho, rthk: the refined layered velocity model
  !   !!         
  !   implicit none
  !   integer NL
  !   parameter (NL=200)
  !   integer mmax,rmax
  !   real minthk0
  !   real minthk
  !   real dep(*),vp(*),vs(*),rho(*)
  !   real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
  !   integer nsublay(NL)
  !   real thk,newthk,initdep 
  !   integer i,j,k,ngrid
  
  !   k = 0
  !   initdep = 0.0
  !   do i = 1, mmax-1
  !     thk = dep(i+1)-dep(i)
  !     minthk = thk/minthk0
  !     nsublay(i) = int((thk+1.0e-4)/minthk) + 1
  !     ngrid = nsublay(i)+1
  !     newthk = thk/nsublay(i)
  !     do j = 1, nsublay(i)
  !       k = k + 1
  !       rthk(k) = newthk
  !       rdep(k) = initdep + rthk(k)
  !       initdep = rdep(k)
  !       rvp(k) = vp(i)+(2*j-1)*(vp(i+1)-vp(i))/(2*nsublay(i))
  !       rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*nsublay(i))
  !       rrho(k) = rho(i)+(2*j-1)*(rho(i+1)-rho(i))/(2*nsublay(i))
  !     enddo
  !   enddo
  !   !! half space model
  !   k = k + 1
  !   rthk(k) = 0.0
  !   rvp(k) = vp(mmax)
  !   rvs(k) = vs(mmax)
  !   rrho(k) = rho(mmax)        
  !   rdep(k) = dep(mmax)
  
  !   rmax = k
  
  !   !!       do i = 1, mmax
  !   !!          write(*,*) dep(i),vp(i),vs(i),rho(i)
  !   !!       enddo
  !   !!       print *, '---------------------------------'
  !   !!       do i = 1, rmax
  !   !!          write(*,*) rdep(i),rthk(i),rvp(i),rvs(i),rrho(i)
  !   !!       enddo
  
  !   return
  !   end