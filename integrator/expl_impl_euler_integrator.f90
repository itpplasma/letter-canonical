module expl_impl_euler_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: integrator_t

    implicit none

    type, extends(integrator_t) :: expl_impl_euler_integrator_t
        real(dp) :: rmu, ro0, rtol
        contains
        procedure :: timestep
    end type expl_impl_euler_integrator_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(expl_impl_euler_integrator_t), intent(in) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        real(dp) :: tstart, tend

        ierr = 0
        tstart = 0d0
        tend = dtau

        call odeint_allroutines(z, 5, tstart, tend, self%rtol, ydot)

        contains

        subroutine ydot(tau, y, vy)
            use velo_sub, only: velo

            real(dp), intent(in) :: tau
            real(dp), dimension(5), intent(in) :: y
            real(dp), dimension(5), intent(inout) :: vy

            call velo(tau, y, vy, self%rmu, self%ro0, magfie_can)

        end subroutine ydot
    end subroutine timestep

end module expl_impl_euler_integrator
