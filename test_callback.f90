program test_callback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use callback, only: cut_callback_t
    implicit none

    type(cut_callback_t) :: cut_callback

    real(dp) :: t_test(9) = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0]
    real(dp) :: z_test(2, 9) = 0d0

    real(dp) :: t_event, z_event(2)
    real(dp) :: t_expected = 4.5d0
    real(dp) :: z_expected(2) = [0d0, 0.5d0]

    integer :: i

    cut_callback%distance => zero_crossed_z1
    cut_callback%event => record_event

    z_test(1,:) = [-3.5d0, -2.5d0, -1.5d0, -0.5d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0]
    z_test(2,:) = [-3d0, -2d0, -1d0, 0d0, 1d0, 2d0, 3d0, 4d0, 5d0]

    do i = 1, 9
        call cut_callback%execute(t_test(i), z_test(:,i))
    end do

    if (abs(t_event - t_expected) > 1d-12) error stop
    if (any(abs(z_event - z_expected) > 1d-12)) error stop

    print *, "test_callback: OK"

    contains

    function zero_crossed_z1(t, z) result(distance)
        real(dp), intent(in) :: t, z(:)
        real(dp) :: distance

        distance = z(1)
    end function zero_crossed_z1


    subroutine record_event(t, z)
        real(dp), intent(in) :: t, z(:)

        t_event = t
        z_event = z
    end subroutine record_event
end program test_callback
