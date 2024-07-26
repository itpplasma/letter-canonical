program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use canonical, only: twopi
    use letter_canonical, only: init, stop_on_error, trace_orbit, write_output

    implicit none

    character(1024) :: input_file = "letter_canonical.in"

    call print_usage

    if (command_argument_count() > 0) then
        call get_command_argument(1, input_file)
    end if

    call init(input_file)
    call stop_on_error

    call trace_orbit
    call stop_on_error

    call write_output
    call stop_on_error

    contains


    subroutine print_usage
        print *, "Usage: letter_canonical.x input_file R0 phi0 Z0 vpar0"
    end subroutine print_usage


end program main
