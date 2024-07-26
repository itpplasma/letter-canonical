program test_orbit
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use omp_lib
  use interpolate, only: evaluate_splines_3d
  use magfie, only: field_t
  use magfie_factory, only: magfie_type_from_string
  use canonical, only: twopi, init_canonical, init_transformation, &
    init_canonical_field_components, init_splines_with_psi, spl_R_of_xc
  use field_can, only: field_can_t, field_can_cyl_t, field_can_albert_t
  use integrator
  use rk45_can_integrator, only: magfie_can
  use odeint_sub

  implicit none

  integer, parameter :: n_r=100, n_phi=64, n_z=75
  integer :: nt
  character(*), parameter :: outname = "orbit.out"

  class(field_t), allocatable :: field_type

  real(dp) :: rmin, rmax, zmin, zmax
  real(dp) :: z0(5), starttime, endtime, tau1, tau2
  real(dp), allocatable :: z(:, :)

  integer :: kt, nmax

  character(16) :: magfie_type="tok"

  real(dp) :: rmu=1d5, ro0=20000d0*1d0  ! 20000 Gauss, 1cm Larmor radius

  real(dp) :: R

  nt = 8000

  ! Workaround, otherwise not initialized without perturbation field
  rmin = 75.d0
  rmax = 264.42281879194627d0
  zmin = -150.d0
  zmax = 147.38193979933115d0

  field_type = magfie_type_from_string(trim(magfie_type))

  print *, "init_magfie ..."
  call field_type%init_magfie()

  print *, "init_canonical ..."
  call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
    [rmax, twopi, zmax], field_type)

  print *, "init_transformation ..."
  call init_transformation

  print *, "init_canonical_field_components ..."
  call init_canonical_field_components

  print *, "init_splines_with_psi ..."
  call init_splines_with_psi

  ! Initial conditions
  z0(1) = 15000000d0 ! psi tok
  !z0(1) = -5.5d0 ! psi test
  z0(2) = -76.5d0 ! z
  z0(3) = -2.0d0*3.1415d0  ! phi
  z0(4) = 1.0d0  ! normalized momentum
  z0(5) = 0.0d0  ! normalized pitch parameter

  allocate(z(5,nt))

  z=0d0

  z(:,1) = z0
  starttime = omp_get_wtime()
  nmax = nt
  do kt = 2, nt
    tau1 = 0d0
    tau2 = 1d0
    z(:,kt) = z(:,kt-1)
    call odeint_allroutines(z(:,kt), 5, tau1, tau2, 1d-8, velo_can)
  end do
  endtime = omp_get_wtime()
  print *, trim(outname), endtime-starttime

  open(unit=20, file=outname, action='write', recl=4096)
  do kt = 1, nmax
    call evaluate_splines_3d(spl_R_of_xc, z(1:3,kt), R)
    write(20,*) z(:,kt), R
  end do
  close(20)
  deallocate(z)

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine velo_can(tau,z,vz)
!
!
!  Computes the components of the 5-D velocity          -  vz
!  for given set of phase space coordinates             -  z.
!
!  Warning !!!  The dimensionless time is used (see below)
!
!  Order of the coordinates is the following:
!  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
!                                 dimension of the general covariant
!                                 space coordinate system
!  z(4) = p                     - momentum  module normalized to
!                                 thermal momentum and sqrt(2);
!  z(5) = alambd                - cosine of the pitch-angle
!
!  Input parameters:
!            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
!                     z      -   see above
!            common:  rmu    -   inverse relativistic temperature
!                     ro0    -   Larmor radius for the reference
!                                magnetic field and temperature:
!                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
!  Output parameters:
!            formal:  vz     -   see above
!
!  Called routines: magfie_can, magfie_vmec, elefie_can, magfie_boozer
!

    implicit none

    integer :: i

    double precision tau,z,vz
    double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
    double precision derphi
    double precision p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
    double precision rmumag,rovsqg,rosqgb,rovbm
    double precision a_phi,a_b,a_c,hstar
    double precision s_hc,hpstar,phidot,blodot,bra

    dimension z(5),vz(5)
    dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
    dimension derphi(3)
    dimension a_phi(3),a_b(3),a_c(3),hstar(3)

    do 1 i=1,3
      x(i)=z(i)
1   continue

! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
!            bmod   - dimensionless magnetic field module: bmod=B/B_ref
!            sqrtg  - Jacobian of space coordinates (square root of
!                     metric tensor
!            bder   - derivatives of logarithm of bmod over space coords
!                     (covariant vector)
!            hcovar - covariant components of the unit vector along
!                     the magnetic field
!            hctrvr - contravariant components of the unit vector along
!                     the magnetic field
!            hcurl  - contravariant components of the curl of this vector


    call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

    ! in elefie: x(i)   - space coords (input, see above)
    !            derphi - derivatives of the dimensionless electric potential
    !                     phihat=e*phi/T over space coords (covar. vector)

    !call elefie_can(x,derphi)
    derphi = 0d0

    p=z(4)
    alambd=z(5)

    p2=p*p
    ovmu=2.d0/rmu
    gamma2=p2*ovmu+1.d0
    gamma=dsqrt(gamma2)
    ppar=p*alambd
    ! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
    vpa=ppar/gamma
    coala=(1.d0-alambd**2)
    ! rmumag - magnetic moment
    rmumag=.5d0*p2*coala/bmod

    rovsqg=ro0/sqrtg
    rosqgb=.5d0*rovsqg/bmod
    rovbm=ro0/bmod

    a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
    a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
    a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
    a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
    a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
    a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg

    s_hc=0.d0
    do i=1,3
      a_c(i)=hcurl(i)*rovbm
      s_hc=s_hc+a_c(i)*hcovar(i)
      hstar(i)=hctrvr(i)+ppar*a_c(i)
    enddo
    hpstar=1.d0+ppar*s_hc

    ! velocities in the coordinate space

    ! phidot - derivative of the dmls el. potential over dmls time
    ! blodot - derivative of the logarith of the mag. field module over dmls time
    phidot=0.d0
    blodot=0.d0
    do i=1,3
      bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
      vz(i)=bra/hpstar
      phidot=phidot+vz(i)*derphi(i)
      blodot=blodot+vz(i)*bder(i)
    enddo

    ! velocities in the phase space

    vz(4)=-0.5d0*gamma*phidot/p
    !      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
    vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
      + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))

  end subroutine velo_can

end program test_orbit
