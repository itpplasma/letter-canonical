module rk45_cyl_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: integrator_t

    implicit none

    type, extends(integrator_t) :: rk45_cyl_integrator_t

        real(dp) :: rmu, ro0, rtol

        contains

        procedure :: timestep => timestep
    end type rk45_cyl_integrator_t

    abstract interface
        subroutine field_p(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ, &
                           dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)
            import :: dp
            real(dp), intent(in) :: r,p,z
            real(dp), intent(inout) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ, &
                dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
            real(dp), intent(inout), optional :: Ar,Ap,Az

        end subroutine field_p
    end interface

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(rk45_cyl_integrator_t), intent(in) :: self
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

            call velo(tau, y, vy, self%rmu, self%ro0, magfie_cyl_tok)

        end subroutine ydot
    end subroutine timestep


    subroutine magfie_cyl_tok(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        use magfie_tok, only: my_field

        real(dp), intent(in)  :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie_cyl(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, my_field)
    end subroutine magfie_cyl_tok


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine magfie_cyl(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, field)
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
    !
    !  Extra parameters:  field            - procedure of magnetic field backend

        real(dp), intent(in)  :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        procedure(field_p) :: field

        real(dp) :: hr,hf,hz

        real(dp) :: br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ

        CALL field(x(1),x(2),x(3),br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

        bmod = dsqrt(br**2 + bf**2 + bz**2)  ! B
        sqrtg = x(1)
        hr = br/bmod
        hf = bf/bmod
        hz = bz/bmod

        bder(1) = (brr*hr + bfr*hf + bzr*hz) / bmod
        bder(2) = (brf*hr + bff*hf + bzf*hz) / bmod
        bder(3) = (brz*hr + bfz*hf + bzz*hz) / bmod

        hcovar(1) = hr
        hcovar(2) = hf*x(1)
        hcovar(3) = hz

        hctrvr(1) = hr
        hctrvr(2) = hf/x(1)
        hctrvr(3) = hz

        ! hcurl =
        hcurl(1)=((bzf-x(1)*bfz)/bmod    + hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
        hcurl(2)=((brz-bzr)/bmod         + hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
        hcurl(3)=((bf+x(1)*bfr-brf)/bmod + hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
        return
    end subroutine magfie_cyl
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



end module rk45_cyl_integrator
