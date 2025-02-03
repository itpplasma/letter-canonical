module rk45_cyl_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: integrator_t
    use magfie_cyl_sub, only: field_p, magfie_cyl_tok, n_field_evaluations

    implicit none


    type, extends(integrator_t) :: rk45_cyl_integrator_t
        real(dp) :: rmu, ro0, rtol
        contains
        procedure :: timestep
        procedure :: get_field_evaluations
    end type rk45_cyl_integrator_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(rk45_cyl_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        real(dp) :: tstart, tend, rmu, ro0

        ierr = 0
        tstart = 0d0
        tend = dtau
        rmu = self%rmu
        ro0 = self%ro0

        call odeint_allroutines(z, 5, tstart, tend, self%rtol, ydot)

        contains

        subroutine ydot(tau, y, vy)
            use velo_sub, only: velo

            real(dp), intent(in) :: tau
            real(dp), dimension(5), intent(in) :: y
            real(dp), dimension(5), intent(inout) :: vy

            call velo(tau, y, vy, rmu, ro0, magfie_cyl_tok)

        end subroutine ydot
    end subroutine timestep

    function get_field_evaluations(self) result(n)
        class(rk45_cyl_integrator_t), intent(in) :: self
        integer(8) :: n

        n = n_field_evaluations
    end function get_field_evaluations


end module rk45_cyl_integrator
