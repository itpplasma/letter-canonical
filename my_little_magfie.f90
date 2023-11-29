module my_little_magfie
    implicit none

contains

    subroutine eval_field_B(x, B)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3)

        real(8) :: dummy

        call field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
            dummy,dummy,dummy,dummy,dummy,dummy)
    end subroutine eval_field_B

    subroutine eval_field_B_and_A(x, y, z, Bx, By, Bz, Ax, Ay, Az)
        real(8), intent(in) :: x, y, z
        real(8), intent(inout) :: Bx, By, Bz, Ax, Ay, Az

        Bx = 0.0
        By = 0.0
        Bz = 2.0

        Ax = 0.0
        Ay = 3.0
        Az = 0.0
    end subroutine eval_field_B_and_A
end module my_little_magfie
