module integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can

    use integrator_base
    use integrator_euler0
    use integrator_euler1
    use integrator_rk45
    use rk45_cyl_integrator, only: rk45_cyl_integrator_t
    use rk45_can_integrator, only: rk45_can_integrator_t
    use expl_impl_euler_integrator, only: expl_impl_euler_integrator_t

    implicit none

    type :: integrator_config_t
        character(len=64) :: integ_type
        character(len=64) :: spatial_coords
        character(len=64) :: momentum_coord
    end type integrator_config_t

    contains

    function create_integrator(config) result(integ)
        class(integrator_t), allocatable :: integ

        type(integrator_config_t), intent(in) :: config

        if (config%momentum_coord == "pphi") then
            integ = create_integrator_can(config)
            return
        end if

    end function create_integrator


    function create_integrator_can(config) result(integ)
        class(integrator_t), allocatable :: integ

        type(integrator_config_t), intent(in) :: config

        class(field_can_t), allocatable :: field
        class(symplectic_integrator_t), allocatable :: sympl_integ

        field = create_field_can(config%spatial_coords)

        select case(config%integ_type)
            case("expl_euler")
                sympl_integ = symplectic_integrator_euler0_t(field)
            case("expl_impl_euler")
                sympl_integ = symplectic_integrator_euler1_t(field)
            case("rk45")
                sympl_integ = symplectic_integrator_rk45_t(field)
            case default
                print *, "create_integrator_can: Unknown type ", config%integ_type
                error stop
        end select
    end function create_integrator_can

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
