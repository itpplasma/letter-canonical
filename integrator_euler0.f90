module integrator_euler0
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use integrator_base, only: symplectic_integrator_t, symplectic_integrator_data_t
    use field_can, only: field_can_t, field_can_data_t, get_derivatives2

    implicit none

    type, extends(symplectic_integrator_t) :: symplectic_integrator_euler0_t
        class(field_can_t), allocatable :: field

        contains

        procedure :: timestep => orbit_timestep_euler0
    end type symplectic_integrator_euler0_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine orbit_timestep_euler0(self, si, f, ierr)
    !
        class(symplectic_integrator_euler0_t), intent(in) :: self
        type(symplectic_integrator_data_t), intent(inout) :: si
        type(field_can_data_t), intent(inout) :: f
        integer, intent(out) :: ierr

        integer, parameter :: n = 2

        integer :: ktau

        ierr = 0
        ktau = 0
        do while(ktau .lt. si%ntau)

            call self%field%evaluate(f, si%z(1), si%z(2), si%z(3), 2)
            call get_derivatives2(f, si%z(4))

            si%z(1) = si%z(1) - si%dt*(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
            si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
            si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph
            si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

            ktau = ktau+1
        enddo
    end subroutine orbit_timestep_euler0



end module integrator_euler0
