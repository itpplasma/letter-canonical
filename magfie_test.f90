module magfie_test
    use magfie, only: FieldType

    implicit none

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

        AR = 0.0d0
        Aphi = 0.5d0*R
        AZ = -log(R)

        BR = 0.0d0
        Bphi = 1.0d0/R
        BZ = 1.0d0
    end subroutine compute_abfield

    subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
        class(TestFieldType), intent(in) :: self
        real(8), intent(in) :: R, phi, Z
        real(8), intent(out) :: BR, Bphi, BZ

        BR = 0.0d0
        Bphi = 1.0d0/R
        BZ = 1.0d0
    end subroutine compute_bfield

end module magfie_test
