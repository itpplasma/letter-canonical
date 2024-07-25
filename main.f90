program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use canonical, only: twopi
    use callback, only: callback_pointer_t, cut_callback_t
    use letter_canonical, only: init, stop_on_error, trace_orbit, write_output, &
        from_internal_coordinates

    implicit none

    character(1024) :: input_file = "letter_canonical.in"
    character(1024) :: arg

    real(dp) :: R0=162.6d0, phi0=-6.283d0, Z0=-76.5d0, vpar0=1d0

    real(dp), allocatable :: z_out(:,:)

    class(callback_pointer_t), allocatable :: callbacks(:)
    class(cut_callback_t), allocatable :: cut_callback

    call print_usage

    if (command_argument_count() > 0) then
        call get_command_argument(1, input_file)
    end if

    if (command_argument_count() > 4) then
        call get_command_argument(2, arg)
        read(arg, *) R0
        call get_command_argument(3, arg)
        read(arg, *) phi0
        call get_command_argument(4, arg)
        read(arg, *) Z0
        call get_command_argument(5, arg)
        read(arg, *) vpar0
    end if

    call init(input_file)
    call stop_on_error

    allocate(callbacks(1))
    allocate(cut_callback)
    cut_callback%distance => zero_crossed_phi
    cut_callback%event => print_state
    allocate(callbacks(1)%item, source=cut_callback)

    call trace_orbit([R0, phi0, Z0, 1d0, vpar0], z_out, callbacks)
    call stop_on_error

    call write_output(z_out)
    call stop_on_error

    contains

    function zero_crossed_phi(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = -modulo(z(2), twopi) + 0.5d0*twopi
    end function zero_crossed_phi


    subroutine print_state(t, z)
        real(dp), intent(in) :: t, z(:)
        write(999, *) t, z
    end subroutine print_state


    subroutine print_usage
        print *, "Usage: letter_canonical.x input_file R0 phi0 Z0 vpar0"
    end subroutine print_usage


end program main
