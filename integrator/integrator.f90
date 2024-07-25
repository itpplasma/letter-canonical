module integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can

    use integrator_base
    use integrator_euler0
    use integrator_euler1
    use integrator_rk45
    use rk45_cyl_integrator
    use rk45_can_integrator
    use expl_impl_euler_integrator_t

    implicit none

    contains

    function create_integrator(integrator_type, field) result(integ)
        character(*), intent(in) :: integrator_type
        class(field_can_t), allocatable :: field
        class(symplectic_integrator_t), allocatable :: integ

        select case(integrator_type)
            case("euler0")
                integ = symplectic_integrator_euler0_t(field)
            case("euler1")
                integ = symplectic_integrator_euler1_t(field)
            case("rk45")
                integ = symplectic_integrator_rk45_t(field)
            case default
                print *, "create_integrator: Unknown integrator type ", integrator_type
                error stop
        end select
    end function create_integrator

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
