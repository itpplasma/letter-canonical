module expl_impl_euler_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use field_can, only: field_can_t, field_can_data_t
    use integrator_base, only: integrator_t, symplectic_integrator_data_t
    use integrator_euler1, only: symplectic_integrator_euler1_t

    implicit none

    type, extends(integrator_t) :: expl_impl_euler_integrator_t
        type(symplectic_integrator_data_t) :: si
        type(field_can_data_t) :: f
        class(symplectic_integrator_euler1_t), allocatable :: sym_integ
        contains
        procedure :: init
        procedure :: timestep
    end type expl_impl_euler_integrator_t

    contains

    subroutine init(self, field)
        class(expl_impl_euler_integrator_t), intent(inout) :: self
        class(field_can_t), intent(in) :: field

        self%sym_integ = symplectic_integrator_euler1_t(field)
    end subroutine init

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(expl_impl_euler_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        self%si%z = z
        call self%sym_integ%timestep(self%si, self%f, ierr)
        z(1:4) = self%si%z
    end subroutine timestep

end module expl_impl_euler_integrator
