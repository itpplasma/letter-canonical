module integrator_rk45
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use integrator_base, only: symplectic_integrator_t, symplectic_integrator_data_t
    use field_can, only: field_can_t, field_can_data_t, get_derivatives2

    implicit none

    type, extends(symplectic_integrator_t) :: symplectic_integrator_rk45_t
        class(field_can_t), allocatable :: field

        contains

        procedure :: timestep => orbit_timestep_rk45
    end type symplectic_integrator_rk45_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine orbit_timestep_rk45(self, si, f, ierr)
    !
        class(symplectic_integrator_rk45_t), intent(in) :: self
        type(symplectic_integrator_data_t), intent(inout) :: si
        type(field_can_data_t), intent(inout) :: f
        integer, intent(out) :: ierr

        integer, parameter :: ndim = 4

        integer :: ktau
        real(dp) :: tstart, tend

        ierr = 0
        ktau = 0
        do while(ktau .lt. si%ntau)
            tstart = ktau*si%dt
            tend = (ktau+1)*si%dt

            call odeint_allroutines(si%z, ndim, tstart, tend, si%rtol, velo)

            ktau = ktau+1
        enddo

        contains

        subroutine velo(t, y, dy)
            use magfie, only: compute_abfield

            real(dp), intent(in) :: t
            real(dp), dimension(4), intent(in) :: y
            real(dp), dimension(4), intent(inout) :: dy

            call self%field%evaluate(f, si%z(1), si%z(2), si%z(3), 2)
            call get_derivatives2(f, si%z(4))

            dy(1) = -(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
            dy(2) = f%dH(1)/f%dpth(1)
            dy(3) = (f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph
            dy(4) = -(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

        end subroutine velo
    end subroutine orbit_timestep_rk45


end module integrator_rk45
