program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use callback, only: callback_t, passing_cut_callback_t
    use letter_canonical, only: init, stop_on_error, trace_orbit, write_output

    implicit none

    character(1024) :: input_file = "letter_canonical.in"
    character(1024) :: arg

    real(dp) :: R0=162.6d0, phi0=-6.283d0, Z0=-76.5d0, vpar0=0d0

    real(dp), allocatable :: z_out(:,:)

    class(callback_t), allocatable :: callbacks(:)

    allocate(passing_cut_callback_t :: callbacks(1))
    allocate(passing_cut_callback_t :: callbacks(2))

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

    call trace_orbit([R0, phi0, Z0, 1d0, vpar0], z_out, callbacks)
    call stop_on_error

    call write_output(z_out)
    call stop_on_error

    contains

    subroutine print_usage
        print *, "Usage: letter_canonical.x input_file R0 phi0 Z0 vpar0"
    end subroutine print_usage

end program main
