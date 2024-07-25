module callback_base
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type, abstract :: callback_t
        contains
        procedure(execute), deferred :: execute
    end type callback_t


    abstract interface
        subroutine execute(self, t, z)
            import callback_t, dp
            class(callback_t), intent(inout) :: self
            real(dp), intent(in) :: t
            real(dp), intent(in) :: z(:)
        end subroutine execute
    end interface


    type :: callback_pointer_t
    ! Pointer to a callback object, required for a list of callbacks
        class(callback_t), allocatable :: item
        contains
        procedure :: execute => execute_item
    end type callback_pointer_t

    contains

    subroutine execute_item(self, t, z)
        class(callback_pointer_t), intent(inout) :: self
        real(dp), intent(in) :: t
        real(dp), intent(in) :: z(:)

        call self%item%execute(t, z)
    end subroutine execute_item

end module callback_base
