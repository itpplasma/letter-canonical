module magfie
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type, abstract :: field_t
        contains
        procedure(init_magfie), deferred :: init_magfie
        procedure(compute_abfield), deferred :: compute_abfield
        procedure(compute_bfield), deferred :: compute_bfield
    end type field_t

    interface
        subroutine init_magfie(self)
            import :: field_t
            class (field_t), intent(in) :: self
        end subroutine
    end interface

    interface
        subroutine compute_abfield(self, R, phi, Z, AR, Aphi, AZ, BR, Bphi, BZ)
            import :: field_t, dp
            class (field_t), intent(in) :: self
            real(dp), intent(in) :: R, phi, Z
            real(dp), intent(out) :: AR, Aphi, AZ, BR, Bphi, BZ
        end subroutine
    end interface

    interface
        subroutine compute_bfield(self, R, phi, Z, BR, Bphi, BZ)
            import :: field_t, dp
            class (field_t), intent(in) :: self
            real(dp), intent(in) :: R, phi, Z
            real(dp), intent(out) :: BR, Bphi, BZ
        end subroutine
    end interface

end module magfie
