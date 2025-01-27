module magfie_cyl_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    integer(8) :: n_field_evaluations = 0

    !$omp threadprivate(n_field_evaluations)

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

    subroutine magfie_cyl_tok(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        use magfie_tok, only: my_field

        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie_cyl(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, my_field)
        n_field_evaluations = n_field_evaluations + 1
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
        !                     bder             - covariant components of (grad B)/B
        !                     hcovar           - covariant components of \bB/B
        !                     hctrvr           - contravariant components of \bB/B
        !                     hcurl            - contravariant components of curl (\bB/B)
        !
        !  Extra parameters:  field            - procedure of magnetic field backend

        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        procedure(field_p) :: field

        real(dp) :: hr, hf, hz

        real(dp) :: br, bf, bz, BRR, BRF, BRZ, BFR, BFF, BFZ, BZR, BZF, BZZ

        call field(x(1), x(2), x(3), br, bf, bz, BRR, BRF, BRZ, BFR, BFF, BFZ, BZR, BZF, BZZ)

        bmod = dsqrt(br**2 + bf**2 + bz**2)  ! B
        sqrtg = x(1)
        hr = br/bmod
        hf = bf/bmod
        hz = bz/bmod

        bder(1) = (brr*hr + bfr*hf + bzr*hz)/bmod
        bder(2) = (brf*hr + bff*hf + bzf*hz)/bmod
        bder(3) = (brz*hr + bfz*hf + bzz*hz)/bmod

        hcovar(1) = hr
        hcovar(2) = hf*x(1)
        hcovar(3) = hz

        hctrvr(1) = hr
        hctrvr(2) = hf/x(1)
        hctrvr(3) = hz

        ! hcurl =
        hcurl(1) = ((bzf - x(1)*bfz)/bmod + hcovar(2)*bder(3) - hcovar(3)*bder(2))/sqrtg
        hcurl(2) = ((brz - bzr)/bmod + hcovar(3)*bder(1) - hcovar(1)*bder(3))/sqrtg
        hcurl(3) = ((bf + x(1)*bfr - brf)/bmod + hcovar(1)*bder(2) - hcovar(2)*bder(1))/sqrtg
        return
    end subroutine magfie_cyl

end module magfie_cyl_sub
