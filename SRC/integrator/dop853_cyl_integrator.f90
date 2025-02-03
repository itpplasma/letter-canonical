module dop853_cyl_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use integrator_base, only: integrator_t
    use magfie_cyl_sub, only: field_p, magfie_cyl_tok, n_field_evaluations
    use dop853_module, only: dop853_class

    implicit none

    type, extends(integrator_t) :: dop853_cyl_integrator_t
        real(dp) :: rmu, ro0, rtol
    contains
        procedure :: timestep
        procedure :: get_field_evaluations
    end type dop853_cyl_integrator_t

    logical :: firstcall = .true.
    type(dop853_class) :: dop853_backend

contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep(self, z, dtau, ierr)
    !
        class(dop853_cyl_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: dtau
        integer, intent(out) :: ierr

        real(dp) :: tstart, tend, rmu, ro0
        real(dp) :: rtol(5), atol(5)
        logical :: status_ok
        integer :: idid

        ierr = 0
        tstart = 0d0
        tend = dtau
        rmu = self % rmu
        ro0 = self % ro0

        if (firstcall) then
            call dop853_backend%initialize(n=5, fcn=ydot, status_ok=status_ok)
            if (.not. status_ok) error stop 'dop853 initialization error'
            firstcall = .false.
        end if

        rtol(:) = self%rtol
        atol(:) = 1d-13
        call dop853_backend%integrate(tstart, z, tend, rtol, atol, 0, idid)

    contains

        subroutine ydot(me, x, y, f)
            use velo_sub, only: velo

            class(dop853_class), intent(inout) :: me
            real(dp), intent(in) :: x
            real(dp), dimension(:), intent(in) :: y
            real(dp), dimension(:), intent(out) :: f

            call velo(x, y, f, rmu, ro0, magfie_cyl_tok)

        end subroutine ydot
    end subroutine timestep

    function get_field_evaluations(self) result(n)
        class(dop853_cyl_integrator_t), intent(in) :: self
        integer(8) :: n

        n = n_field_evaluations
    end function get_field_evaluations

end module dop853_cyl_integrator
