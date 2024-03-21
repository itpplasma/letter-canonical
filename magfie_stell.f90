module magfie
    implicit none

    contains

    subroutine init_magfie
    end subroutine init_magfie

    subroutine compute_bfield(R, phi, Z, BR, Bphi, BZ)
        real(8), intent(in) :: R, phi, Z
        real(8), intent(out) :: BR, Bphi, BZ

        BR = 0.0d0
        Bphi = 1.0d0/R
        BZ = 0.0d0
    end subroutine compute_bfield

    subroutine compute_abfield(R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
        real(8), intent(in) :: R, phi, Z
        real(8), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ

        AR = 1.0d0
        Aphi = 1.0d0
        AZ = -log(R)

        BR = 0.0d0
        Bphi = 1.0d0/R
        BZ = 0.0d0
    end subroutine compute_abfield

end module magfie
