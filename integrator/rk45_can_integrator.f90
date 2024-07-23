module rk45_can_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: integrator_t

    implicit none

    type, extends(integrator_t) :: rk45_can_integrator_t

        real(dp) :: rmu, ro0, rtol

        contains

        procedure :: timestep => timestep
    end type rk45_can_integrator_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(rk45_can_integrator_t), intent(in) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        real(dp) :: tstart, tend

        ierr = 0
        tstart = 0d0
        tend = dtau

        call odeint_allroutines(z, 5, tstart, tend, self%rtol, ydot)

        contains

        subroutine ydot(tau, y, vy)
            use velo_sub, only: velo

            real(dp), intent(in) :: tau
            real(dp), dimension(5), intent(in) :: y
            real(dp), dimension(5), intent(inout) :: vy

            call velo(tau, y, vy, self%rmu, self%ro0, magfie_can)

        end subroutine ydot
    end subroutine timestep


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    !  Computes magnetic field and derivatives with bmod in units of the magnetic code
    !
    !  Input parameters:
    !            formal:  x                - array of cylindrical coordinates R, phi, Z
    !  Output parameters:
    !            formal:  bmod             - magnetic field module
    !                     sqrtg            - metric determinant
    !                     bder             - covariant components of (grad B)/Bmod
    !                     hcovar           - covariant components of B/Bmod
    !                     hctrvr           - contravariant components of B/Bmod
    !                     hcurl            - contravariant components of (curl B)/Bmod

        use interpolate, only: evaluate_splines_3d_der2
        use canonical, only: spl_Aphi_of_xc, spl_hph_of_xc, spl_hth_of_xc, spl_Bmod_of_xc

        real(dp), intent(in) :: x(3)

        real(dp), intent(out) :: bmod, sqrtg
        real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

        integer, parameter :: reorder(3) = [1, 3, 2]  ! r, ph, th -> r, th, ph

        real(dp) :: aph, hph, hth, sqrtg_bmod
        real(dp), dimension(3) :: daph, dhph, dhth, dbmod
        real(dp), dimension(6) :: dummy

        call evaluate_splines_3d_der2(spl_Aphi_of_xc, x(reorder), aph, daph, dummy)
        call evaluate_splines_3d_der2(spl_Bmod_of_xc, x(reorder), bmod, dbmod, dummy)
        call evaluate_splines_3d_der2(spl_hth_of_xc, x(reorder), hth, dhth, dummy)
        call evaluate_splines_3d_der2(spl_hph_of_xc, x(reorder), hph, dhph, dummy)

        daph = daph(reorder)
        dbmod = dbmod(reorder)
        dhth = dhth(reorder)
        dhph = dhph(reorder)

        sqrtg_bmod = hph - hth*daph(1)
        sqrtg = sqrtg_bmod/bmod
        bder = dbmod/bmod

        hcovar(1) = 0.d0
        hcovar(2) = hth
        hcovar(3) = hph

        hctrvr(1) = daph(2)/sqrtg_bmod
        hctrvr(2) = -daph(1)/sqrtg_bmod
        hctrvr(3) = 1.d0/sqrtg_bmod

        hcurl(1) = (dhph(2) - dhth(3))/sqrtg
        hcurl(2) = -dhph(1)/sqrtg
        hcurl(3) = dhth(1)/sqrtg

    end subroutine magfie_can



end module rk45_can_integrator
