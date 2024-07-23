module rk45_cyl_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: integrator_t

    implicit none

    type, extends(integrator_t) :: rk45_cyl_integrator_t

        contains

        procedure :: timestep => timestep
    end type rk45_cyl_integrator_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(rk45_cyl_integrator_t), intent(in) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        integer, parameter :: ndim = 4

        integer :: ktau
        real(dp) :: tstart, tend

        ierr = 0
        ktau = 0

        contains

        subroutine velo(t, y, dy)

            real(dp), intent(in) :: t
            real(dp), dimension(4), intent(in) :: y
            real(dp), dimension(4), intent(inout) :: dy


        end subroutine velo
    end subroutine timestep


end module rk45_cyl_integrator
