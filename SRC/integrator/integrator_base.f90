module integrator_base
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_data_t

    implicit none

    type, abstract :: integrator_t
        contains
        procedure(timestep), deferred :: timestep
        procedure(get_field_evaluations), deferred :: get_field_evaluations
    end type integrator_t

    abstract interface
        subroutine timestep(self, z, dtau, ierr)
            import integrator_t, dp
            class(integrator_t), intent(inout) :: self
            real(dp), intent(inout) :: z(:)
            real(dp), intent(in) :: dtau
            integer, intent(out) :: ierr
        end subroutine timestep
    end interface

    abstract interface
        function get_field_evaluations(self)
            import integrator_t, field_can_data_t, dp
            integer(8) :: get_field_evaluations
            class(integrator_t), intent(in) :: self
        end function get_field_evaluations
    end interface


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
        procedure(symplectic_timestep), deferred :: timestep
    end type symplectic_integrator_t

    abstract interface
        subroutine symplectic_timestep(self, si, f, ierr)
            import symplectic_integrator_t, symplectic_integrator_data_t, &
                field_can_data_t
            class(symplectic_integrator_t), intent(in) :: self
            type(symplectic_integrator_data_t), intent(inout) :: si
            type(field_can_data_t), intent(inout) :: f
            integer, intent(out) :: ierr
        end subroutine symplectic_timestep
    end interface

    contains

    subroutine odeint_rk4(y, nvar, x1, x2, derivs)
        integer, intent(in) :: nvar
        real(dp), intent(in) :: x1, x2
        real(dp), dimension(nvar), intent(inout) :: y
        external :: derivs

        real(dp) :: h, x, hh, h6
        real(dp), dimension(nvar) :: dydx, dyt, dym, yt
        integer :: i

        h = (x2 - x1) / 2d0
        hh = h * 0.5d0
        h6 = h / 6.0d0
        x = x1

        call derivs(x, y, dydx)

        do i = 1, nvar
          yt(i) = y(i) + hh * dydx(i)
        end do
        call derivs(x + hh, yt, dyt)

        do i = 1, nvar
          yt(i) = y(i) + hh * dyt(i)
        end do
        call derivs(x + hh, yt, dym)

        do i = 1, nvar
          yt(i) = y(i) + h * dym(i)
          dym(i) = dym(i) + dyt(i)
        end do
        call derivs(x + h, yt, dyt)

        do i = 1, nvar
          y(i) = y(i) + h6 * (dydx(i) + dyt(i) + 2.0d0 * dym(i))
        end do

      end subroutine odeint_rk4

end module integrator_base
