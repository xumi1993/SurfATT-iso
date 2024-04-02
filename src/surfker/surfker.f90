module surfker

  use constants
  use RayleighWaveKernel
  use utils
  implicit none
  integer, parameter :: iflsph=1, mode=1
  ! integer, parameter :: dp = 8
contains

  subroutine depthkernel_mpi(vel3d,igrid,istart,iend,iwave,igr,tRc,depz,&
                            sen_vsRc,sen_vpRc,sen_rhoRc)
    real(kind=dp), dimension(:,:,:), intent(in) :: vel3d
    real(kind=dp), dimension(:), intent(in) :: tRc
    real(kind=dp), dimension(:), intent(in) :: depz
    real(kind=dp), dimension(:,:,:,:), intent(inout) :: sen_vsRc, sen_vpRc, sen_rhoRc
    integer, intent(in) :: iwave, igr, istart,iend
    integer, dimension(:,:), intent(in) :: igrid
    real(kind=dp), dimension(:), allocatable :: vsz
    real, dimension(:), allocatable :: rvp, rvs, rrho, rthk
    integer :: nx, ny, nz,kmaxRc, ii

    nx = size(vel3d,1)
    ny = size(vel3d,2)
    nz = size(vel3d,3)
    kmaxRc = size(tRc)
    allocate(vsz(nz))
    do ii = istart, iend
      vsz(1:nz)=vel3d(igrid(ii,1),igrid(ii,2),1:nz)
      call depthkernel1d(vsz,nz,iwave,igr,kmaxRc,tRc,depz,&
                    sen_vsRc(1:kmaxRc,igrid(ii,1),igrid(ii,2),1:nz),&
                    sen_vpRc(1:kmaxRc,igrid(ii,1),igrid(ii,2),1:nz),&
                    sen_rhoRc(1:kmaxRc,igrid(ii,1),igrid(ii,2),1:nz))
    enddo
  end subroutine depthkernel_mpi

  subroutine depthkernel1d(vel,nz,iwave,igr,kmaxRc,tRc,depz,&
                          sen_vsRc,sen_vpRc,sen_rhoRc)
    integer ::  nz, kmaxRc, iwave,igr
    real(kind=dp), dimension(nz), intent(in) :: vel, depz
    real(kind=dp), dimension(kmaxRc), intent(in) :: tRc
    real(kind=dp), dimension(kmaxRc,nz), intent(inout) :: sen_vsRc, sen_vpRc, sen_rhoRc
    real(kind=dp), dimension(kmaxRc) :: cp,cg
    real, dimension(nz) :: vpz,vsz,rhoz
    real, dimension(:), allocatable :: rvp,rvs,rrho,rthk
    real(kind=dp), dimension(:), allocatable :: dcdar,dcdbr,dcdhr,dcdrr,ur,uz,tr,tz
    integer :: mmax, i

    vsz = real(vel)
    mmax = nz+1
    allocate(dcdar(mmax), dcdbr(mmax), dcdhr(mmax), dcdrr(mmax))
    dcdar = 0._dp
    dcdbr = 0._dp
    dcdhr = 0._dp
    dcdrr = 0._dp
    allocate(ur(mmax), uz(mmax), tr(mmax), tz(mmax))
    allocate(rvp(mmax), rvs(mmax), rrho(mmax), rthk(mmax))
    call get_vprho(vsz, nz, vpz, rhoz)
    call refinegrid(vpz, vsz, rhoz, real(depz), nz, rvp, rvs, rrho, rthk)
    if (iwave == 2 .and. igr == 0) then
      call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,igr,kmaxRc,&
                      tRc,cp)
      do i = 1, kmaxRc
        call sregn96(rthk,rvp,rvs,rrho,mmax,tRc(i),cp(i),cg(i),&
                  ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph)
                  ! thk,vp,vs,rhom,nlayer,&
                  !   t,cp,cg,dispu,dispw,stressu,stressw,&
                  !   dc2da,dc2db,dc2dh,dc2dr,iflsph)
        sen_vsRc(i,1:nz) = dcdbr(1:nz)
        sen_vpRc(i,1:nz) = dcdar(1:nz)
        sen_rhoRc(i,1:nz) = dcdrr(1:nz)
      enddo
    elseif (iwave == 2 .and. igr == 1) then
      block
        real(kind=dp), dimension(kmaxRc) :: t1, t2, cp, c1, c2, cg
        real(kind=dp), dimension(mmax) :: dudar, dudbr, dudhr, dudrr
        real(kind=dp) :: dt = 0.01
        t1 = tRc * (1+0.5*dt)
        t2 = tRc * (1-0.5*dt)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          tRc,cp)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          t1,c1)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          t2,c2)
        do i = 1, kmaxRc
          call sregnpu(rthk,rvp,rvs,rrho,mmax,tRc(i),cp(i),cg(i),&
                      ur,uz,tr,tz,t1(i),c1(i),t2(i),c2(i),&
                      dcdar,dcdbr,dcdhr,dcdrr,&
                      dudar,dudbr,dudhr,dudrr,iflsph)
                      !  (thk,vp,vs,rhom,nlayer,&
                      !   t,cp,cg,dispu,dispw,stressu,stressw,&
                      !   t1,cp1,t2,cp2,&
                      !   dc2da,dc2db,dc2dh,dc2dr,&
                      !   du2da,du2db,du2dh,du2dr,&
                      !   iflsph)
          sen_vsRc(i,1:nz) = dudbr(1:nz)
          sen_vpRc(i,1:nz) = dudar(1:nz)
          sen_rhoRc(i,1:nz) = dudrr(1:nz)  
        enddo
      end block
    else
      stop 'kernel1D: Only rayleigh wave is supported for now.'
    endif

  end subroutine depthkernel1d

  subroutine fwdsurf3d_mpi(vel, igrid,istart,iend,iwave,igr,tRc,depz,svel)
  !subroutine fwdsurf3d_mpi(vel,nx,ny,nz,igrid,istart,iend,&
  !  iwave,igr,kmaxRc,tRc,depz,minthk,svel)
    real(kind=dp), dimension(:,:,:), intent(in) :: vel
    real(kind=dp), dimension(:), intent(in) :: tRc
    real(kind=dp), dimension(:), intent(in) :: depz
    real(kind=dp), dimension(:,:,:), intent(out) :: svel
    integer, intent(in) :: iwave, igr, istart,iend
    integer, dimension(:,:), intent(in) :: igrid
    real(kind=dp), dimension(:), allocatable :: vsz, tmpv
    integer :: nx, ny, nz,kmaxRc, ii

    nz = size(vel,3)
    kmaxRc = size(tRc)
    allocate(vsz(nz), tmpv(kmaxRc))
    do ii = istart, iend
      vsz(1:nz)=vel(igrid(ii,1),igrid(ii,2),1:nz)
      tmpv = 0._dp
      call fwdsurf1d(vsz,nz,iwave,igr,kmaxRc,tRc,depz,tmpv)
      svel(1:kmaxRc, igrid(ii,1),igrid(ii,2)) = tmpv(1:kmaxRc)
    enddo
  end subroutine fwdsurf3d_mpi

  subroutine fwdsurf1d(vel,nz,iwave,igr,kmaxRc,tRc,depz,svel)
    integer, intent(in) :: nz, kmaxRc, iwave,igr
    real(kind=dp), dimension(nz) :: vel, depz
    real(kind=dp), dimension(kmaxRc) :: tRc
    real(kind=dp), dimension(kmaxRc), intent(out) :: svel

    real, dimension(nz) :: vpz,vsz,rhoz,rthk
    integer mmax
    integer kk

    mmax=nz
    vsz(1:nz)=real(vel(1:nz))
    ! some other emperical relationship maybe better, 
    ! This is from Tomas M.Brocher 2005 BSSA
    call get_vprho(vsz, nz, vpz, rhoz)
    do kk=1,mmax
      rthk(kk) = depz(kk+1)-depz(kk)
    enddo
    !!half space
    rthk(mmax) = 0.
    ! call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax,rdep,&
    !   rvp,rvs,rrho,rthk)
    call surfdisp96(rthk,vpz,vsz,rhoz,mmax,iflsph,iwave,mode,igr,kmaxRc,&
                    tRc,svel)
    ! svel(1:kmaxRc)=cgRc(1:kmaxRc)

  end subroutine fwdsurf1d

  subroutine correct_depth(ker_vs, ker_vp, ker_rho, dip_angle, igrid,istart,iend, depz)
    real(kind=dp), dimension(:,:,:,:), intent(inout) :: ker_vs, ker_vp, ker_rho
    real(kind=dp), dimension(:,:,:), intent(in) :: dip_angle
    integer, intent(in) :: istart, iend
    integer, dimension(:,:), intent(in) :: igrid
    real(kind=dp), dimension(:), intent(in) :: depz
    real(kind=dp), dimension(:), allocatable :: newz
    integer :: i, np, j

    np = size(dip_angle,1)
    do i = istart, iend
      do j = 1, np
        newz = depz/cosd(dip_angle(j, igrid(i,1),igrid(i,2)))
        ker_vs(j,igrid(i,1),igrid(i,2),:) = interp1(depz, ker_vs(j,igrid(i,1),igrid(i,2),:), newz)
        ker_vp(j,igrid(i,1),igrid(i,2),:) = interp1(depz, ker_vp(j,igrid(i,1),igrid(i,2),:), newz)
        ker_rho(j,igrid(i,1),igrid(i,2),:) = interp1(depz, ker_rho(j,igrid(i,1),igrid(i,2),:), newz)
      enddo
    enddo
    
  end subroutine correct_depth

  subroutine refinegrid(vp, vs, rho, dep, nz, rvp, rvs, rrho, rthk)
    integer :: nz
    real, dimension(nz), intent(in) :: vp, vs, rho, dep
    real, dimension(nz+1), intent(out) :: rvp, rvs, rrho, rthk
    integer kk, mmax

    mmax = nz+1
    do kk=1,nz
      rvs(kk) = vs(kk)
      rvp(kk) = vp(kk)
      if (kk == nz) then
        rthk(kk) = rthk(kk-1)
      else
        rthk(kk) = dep(kk+1)-dep(kk)
      endif
      rrho(kk) = rho(kk)
    enddo
    rthk(mmax) = 0.
    rvp(mmax) = rvp(nz)
    rvs(mmax) = rvs(nz)
    rrho(mmax) = rrho(nz)

  end subroutine refinegrid

  subroutine get_vprho(vsz, nz, vpz, rhoz)
    integer :: nz
    real, dimension(nz), intent(in) :: vsz
    real, dimension(nz), intent(out) :: vpz, rhoz
    
    vpz=0.9409 + 2.0947*vsz - 0.8206*vsz**2+ &
        0.2683*vsz**3 - 0.0251*vsz**4

    rhoz=1.6612*vpz - 0.4721*vpz**2 + &
        0.0671*vpz**3 - 0.0043*vpz**4 + & 
        0.000106*vpz**5

  end subroutine 
end module surfker
