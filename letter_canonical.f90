module letter_canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use integrator, only: integrator_t
    implicit none

    integer :: error_code_ = 0

    class(integrator_t), allocatable :: integ

    abstract interface
        subroutine callback_type(z)
            import :: dp
            real(dp), intent(in) :: z(5)
        end subroutine callback_type
    end interface

contains

    subroutine init()
    end subroutine init

    subroutine trace_orbit(z0, dtau, ntau, nskip, z_out, callback)
        real(dp), intent(in) :: z0(5)
        real(dp), intent(in) :: dtau
        integer, intent(in) :: ntau, nskip
        real(dp), intent(inout) :: z_out(:, :)
        procedure(callback_type), optional :: callback

        integer :: kt, nt, ierr
        real(dp) :: z(5)

        nt = ntau/nskip

        if (.not. size(z_out, 2) == nt) then
            call throw_error("trace_orbit: z_out has wrong size")
            return
        end if

        z_out(:, 1) = z0
        do kt = 2, nt
            call integ%timestep(z, dtau, ierr)
            if (ierr /= 0) then
                call throw_error("trace_orbit: error in timestep", ierr)
                return
            end if
            call callback(z)
            if (mod(kt, nskip) == 0) then
                z_out(:, kt/nskip) = z
            end if
        end do
    end subroutine trace_orbit


    subroutine throw_error(msg, error_code)
        character(*), intent(in) :: msg
        integer, intent(in), optional :: error_code
        print *, msg
        error_code_ = error_code
    end subroutine throw_error

end module letter_canonical
