program test_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use magfie, only: FieldType
    use magfie_factory, only: magfie_type_from_string
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components, init_splines_with_psi
    use field_can, only: field_can_t, field_can_cyl_t, field_can_new_t
    use integrator

    implicit none

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer :: nt
    character(*), parameter :: outname = "orbit.out"
    real(dp), parameter :: qe = 1d5, m = 1d0, c = 1d0, mu = 0d0 !1d-5

    class(FieldType), allocatable :: field_type
    class(field_can_t), allocatable :: field
    type(field_can_data_t) :: f
    class(symplectic_integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si

    real(dp) :: rmin, rmax, zmin, zmax
    real(dp) :: z0(4), vpar0, starttime, endtime
    real(dp), allocatable :: out(:, :)

    integer :: kt, ierr, nmax

    ! Configuration in letter_canonical.in
    character(16) :: magfie_type, integrator_type
    namelist /letter_canonical/ magfie_type, integrator_type

    nt = 8000

    ! Workaround, otherwise not initialized without perturbation field
    rmin = 75.d0
    rmax = 264.42281879194627d0
    zmin = -150.d0
    zmax = 147.38193979933115d0

    open(unit=10, file='letter_canonical.in', status='old', action='read')
    read(10, nml=letter_canonical)
    close(10)
    field_type = magfie_type_from_string(trim(magfie_type))

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
    !z0(1) = 173d0   ! r
    z0(2) = 27.5d0 ! z
    z0(3) = -2.0d0*3.1415d0  ! phi
    vpar0 = 0.8d0  ! parallel velocity

    field = field_can_new_t()
    call field_can_init(f, mu, c*m/qe, vpar0)
    call field%evaluate(f, z0(1), z0(2), z0(3), 2)

    z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi
    call get_derivatives2(f, z0(4))
    print *, 'f%hph = ', f%hph
    print *, 'f%Aph = ', f%Aph
    print *, 'f%ro0 = ', f%ro0
    print *, 'vpar0 = ', vpar0
    print *, 'z0 = ', z0

    integ = create_integrator(trim(integrator_type), field)
    call integrator_init(si, field, f, z0, dt=1.0d-2, ntau=1, rtol=1d-8)

    allocate(out(5,nt))

    out=0d0

    out(1:4,1) = z0
    out(5,1) = f%H
    starttime = omp_get_wtime()
    nmax = nt
    do kt = 2, nt
        ierr = 0
        call integ%timestep(si, f, ierr)
        if (.not. ierr==0) then
            print *, si%z
            nmax = kt-1
            print *, 'nmax = ', nmax
            exit
        endif
        out(1:4,kt) = si%z
        out(5,kt) = f%H
    end do
    endtime = omp_get_wtime()
    print *, trim(outname), endtime-starttime

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nmax
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  !
  !
  ! Computes magnetic field module in units of the magnetic code  - bmod,
  ! square root of determinant of the metric tensor               - sqrtg,
  ! derivatives of the logarythm of the magnetic field module
  ! over coordinates                                              - bder,
  ! covariant componets of the unit vector of the magnetic
  ! field direction                                               - hcovar,
  ! contravariant components of this vector                       - hctrvr,
  ! contravariant component of the curl of this vector            - hcurl
  ! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
  ! x(2)=vartheta_c (canonical poloidal angle), x(3)=varphi_c (canonical toroidal angle).
  !
  !  Input parameters:
  !            formal:  x(3)             - array of canonical coordinates
  !  Output parameters:
  !            formal:  bmod
  !                     sqrtg
  !                     bder(3)          - derivatives of $\log(B)$
  !                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
  !                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
  !                     hcurl(3)         - contra-variant components of curl of $\bh$
  !
  !  Called routines: canonical_field
  !
  !
    use canonical, only: spl_Aphi_of_xc, spl_hph_of_xc, spl_hth_of_xc, spl_Bmod_of_xc

    logical :: fullset
    integer :: mode_secders
    double precision :: bmod,sqrtg
    double precision :: r,vartheta_c,varphi_c,                                           &
                        dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                               &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp
    double precision :: Bctr_vartheta,Bctr_varphi,bmod2
    double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl

    integer, parameter :: reorder(3) = [1, 3, 2]  ! dr, dph, dth -> dr, dth, dph
    integer, parameter :: reorder2(6) = [1, 3, 2, 6, 5, 4]
    ! drdr, drdph, drdth, dphdph, dphdth, dthdth ->
    ! drdr, drdth, drdph, dthdth, dthdph, dphdph

    real(dp) :: a, da(3), d2a(6)
  !
    r=x(1)
    vartheta_c=x(2)
    varphi_c=x(3)
  !
    fullset=.false.
    mode_secders=0

    call evaluate_splines_3d_der2(spl_Aphi_of_xc, x(reorder), a, da, d2a)

    dA_phi_dr = a(1)
    d2A_phi_dr2 = d2a(1)

    dA_theta_dr = 1.0d0

    call evaluate_splines_3d_der2(spl_Bmod_of_xc, x(reorder), bmod, da, d2a)

    call evaluate_splines_3d_der2(spl_hph_of_xc, x(reorder), a, da, d2a)


    call evaluate_splines_3d_der2(spl_hth_of_xc, x(reorder), a, da, d2a)

    B_vartheta_c = a(2)*bmod
    d

  !
    sqg_c=r  ! TODO: actual sqrtg, not just rough estimate
    sqrtg=sqg_c
  !

    Bctr_vartheta=-dA_phi_dr/sqg_c
    Bctr_varphi=dA_theta_dr/sqg_c
  !
    bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
    bmod=sqrt(bmod2)
  !
    bder(1)=0.5d0*((dA_theta_dr*dB_varphi_c_dr-dA_phi_dr*dB_vartheta_c_dr-d2A_phi_dr2*B_vartheta_c) &
           /bmod2-dsqg_c_dr)/sqg_c
    bder(2)=0.5d0*((dA_theta_dr*dB_varphi_c_dt-dA_phi_dr*dB_vartheta_c_dt)/bmod2-dsqg_c_dt)/sqg_c
    bder(3)=0.5d0*((dA_theta_dr*dB_varphi_c_dp-dA_phi_dr*dB_vartheta_c_dp)/bmod2-dsqg_c_dp)/sqg_c
  !
    hcovar(1)=0.d0
    hcovar(2)=B_vartheta_c/bmod
    hcovar(3)=B_varphi_c/bmod
  !
    hctrvr(1)=0.d0
    hctrvr(2)=Bctr_vartheta/bmod
    hctrvr(3)=Bctr_varphi/bmod
  !
    hcurl(1)=((dB_varphi_c_dt-dB_vartheta_c_dp)/bmod-bder(2)*hcovar(3)+bder(3)*hcovar(2))/sqg_c
    hcurl(2)=(-dB_varphi_c_dr/bmod+bder(1)*hcovar(3))/sqg_c
    hcurl(3)=(dB_vartheta_c_dr/bmod-bder(1)*hcovar(2))/sqg_c
  !
    end subroutine magfie_can

end program test_orbit
