module callback_base
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type :: callback_t
        procedure(execute), pointer, nopass :: execute => null()
    end type callback_t


    interface
        subroutine execute(self, tau, z)
            import callback_t, dp
            class(callback_t), intent(in) :: self
            real(dp), intent(in) :: tau
            real(dp), intent(inout) :: z(:)
        end subroutine execute
    end interface

end module callback_base
