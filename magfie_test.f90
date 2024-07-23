module magfie_test
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: field_t

    implicit none

    real(dp), parameter :: AMPL = 1.0d-7, AMPL2 = 2.0d-6

    type, extends(field_t) :: Testfield_t
        contains
            procedure :: init_magfie
            procedure :: compute_bfield
            procedure :: compute_abfield
    end type Testfield_t

    contains

    subroutine init_magfie(self)
        class(Testfield_t), intent(in) :: self
    end subroutine init_magfie

    subroutine compute_abfield(self, R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
        class(Testfield_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ

        AR = AMPL*R*cos(phi)
        Aphi = -0.5d0*AMPL2*Z*R
        AZ = -log(R)

        call self%compute_bfield(R, phi, Z, BR, Bphi, BZ)
    end subroutine compute_abfield

    subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
        class(Testfield_t), intent(in) :: self
        real(dp), intent(in) :: R, phi, Z
        real(dp), intent(out) :: BR, Bphi, BZ

        BR = 0.5d0*AMPL2*R
        Bphi = 1.0d0/R
        BZ = -AMPL2*Z + AMPL*sin(phi)
    end subroutine compute_bfield

end module magfie_test
