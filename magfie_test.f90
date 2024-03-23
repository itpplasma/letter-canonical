module magfie_test
    use magfie, only: FieldType

    implicit none

    real(8), parameter :: AMPL = 1.0d-8, AMPL2 = 2.0d-6

    type, extends(FieldType) :: TestFieldType
        contains
            procedure :: init_magfie
            procedure :: compute_bfield
            procedure :: compute_abfield
    end type TestFieldType

    contains

    subroutine init_magfie(self)
        class(TestFieldType), intent(in) :: self
    end subroutine init_magfie

    subroutine compute_abfield(self, R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
        class(TestFieldType), intent(in) :: self
        real(8), intent(in) :: R, phi, Z
        real(8), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ

        AR = AMPL*R*cos(phi)
        Aphi = -0.5d0*AMPL2*Z*R
        AZ = -log(R)

        call self%compute_bfield(R, phi, Z, BR, Bphi, BZ)
    end subroutine compute_abfield

    subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
        class(TestFieldType), intent(in) :: self
        real(8), intent(in) :: R, phi, Z
        real(8), intent(out) :: BR, Bphi, BZ

        BR = 0.5d0*AMPL2*R
        Bphi = 1.0d0/R
        BZ = -AMPL2*Z + AMPL*sin(phi)
    end subroutine compute_bfield

end module magfie_test
