module integrator_base
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_data_t

    implicit none

    type :: symplectic_integrator_data_t
        real(dp) :: atol
        real(dp) :: rtol

        ! Current phase-space coordinates z and old pth
        real(dp), dimension(4) :: z  ! z = (r, th, ph, pphi)
        real(dp) :: pthold

        ! Timestep and variables from z0
        integer :: ntau
        real(dp) :: dt
        real(dp) :: pabs
    end type symplectic_integrator_data_t

    type, abstract :: symplectic_integrator_t
    contains
        procedure(timestep), deferred :: timestep
    end type symplectic_integrator_t

    abstract interface
        subroutine timestep(self, si, f, ierr)
            import symplectic_integrator_t, symplectic_integrator_data_t, &
                field_can_data_t
            class(symplectic_integrator_t), intent(in) :: self
            type(symplectic_integrator_data_t), intent(inout) :: si
            type(field_can_data_t), intent(inout) :: f
            integer, intent(out) :: ierr
        end subroutine timestep
    end interface

end module integrator_base
