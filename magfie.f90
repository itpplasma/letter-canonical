module magfie
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type, abstract :: FieldType
        contains
        procedure(init_magfie), deferred :: init_magfie
        procedure(compute_abfield), deferred :: compute_abfield
        procedure(compute_bfield), deferred :: compute_bfield
    end type FieldType

    interface
        subroutine init_magfie(self)
            import :: FieldType
            class (FieldType), intent(in) :: self
        end subroutine
    end interface

    interface
        subroutine compute_abfield(self, R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
            import :: FieldType, dp
            class (FieldType), intent(in) :: self
            real(dp), intent(in) :: R, phi, Z
            real(dp), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ
        end subroutine
    end interface

    interface
        subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
            import :: FieldType, dp
            class (FieldType), intent(in) :: self
            real(dp), intent(in) :: R, phi, Z
            real(dp), intent(out) :: BR, Bphi, BZ
        end subroutine
    end interface

end module magfie
