module integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can

    use integrator_base
    use integrator_euler1

    implicit none

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine integrator_init(si, field, f, z, dt, ntau, rtol)

        type(symplectic_integrator_data_t), intent(inout) :: si
        class(field_can_t), intent(in) :: field
        type(field_can_data_t), intent(inout) :: f
        double precision, intent(in) :: z(:)
        double precision, intent(in) :: dt
        integer, intent(in) :: ntau
        double precision, intent(in) :: rtol

        si%atol = 1d-15
        si%rtol = rtol

        si%ntau = ntau
        si%dt = dt

        si%z = z

        call field%evaluate(f, z(1), z(2), z(3), 0)
        call get_val(f, si%z(4)) ! for pth
        si%pthold = f%pth
    end subroutine integrator_init


end module integrator
