program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use letter_canonical, only: init, stop_on_error, R0, phi0, Z0, vpar0

    implicit none

    if (command_argument_count() < 4) then
        call print_usage
        stop
    end if

    call init
    call stop_on_error

    contains

    subroutine print_usage
        print *, "Usage: letter_canonical.x input_file R0 phi0 Z0 vpar0"
    end subroutine print_usage

end program main
