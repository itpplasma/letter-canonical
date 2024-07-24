module passing_cut_callback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use callback_base, only: callback_t
    implicit none

    type, extends(callback_t) :: passing_cut_callback_t
        contains
        procedure :: execute => execute
    end type passing_cut_callback_t

    contains

    subroutine execute(self, tau, z)
        class(passing_cut_callback_t), intent(in) :: self
        real(dp), intent(in) :: tau
        real(dp), intent(inout) :: z(:)

        if (z(1) < 0d0) then
            print *, "Cutting at R < 0"
        end if
    end subroutine execute
end module passing_cut_callback
