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
        real(dp) :: zstart(5)
        real(dp) :: dt
        real(dp) :: ro0
        real(dp) :: rtol
        integer :: nskip
    end type integrator_config_t

    contains

    function create_integrator(config, si, f) result(integ)
        class(integrator_t), allocatable :: integ

        type(integrator_config_t), intent(in) :: config
        type(symplectic_integrator_data_t), intent(inout), optional :: si
        type(field_can_data_t), intent(inout), optional :: f

        if (config%momentum_coord == "pphi") then
            integ = create_integrator_pphi(config, si, f)
        else if (config%momentum_coord == "vpar") then
            integ = create_integrator_vpar(config)
        end if

    end function create_integrator


    function create_integrator_pphi(config, si, f) result(integ)
        class(integrator_t), allocatable :: integ

        type(integrator_config_t), intent(in) :: config
        type(symplectic_integrator_data_t), intent(inout) :: si
        type(field_can_data_t), intent(inout) :: f

        class(field_can_t), allocatable :: field
        class(expl_impl_euler_integrator_t), allocatable :: expl_impl_euler_integ

        field = create_field_can(config%spatial_coords)
        call integrator_init(si, field, f, config)

        select case(config%integ_type)
            case("expl_impl_euler")
                expl_impl_euler_integ = expl_impl_euler_integrator_t(si, f)
                call expl_impl_euler_integ%init(field)
                integ = expl_impl_euler_integ
            case default
                print *, "create_integrator_pphi: Unknown type ", config%integ_type
                error stop
        end select
    end function create_integrator_pphi


    function create_integrator_vpar(config) result(integ)
        class(integrator_t), allocatable :: integ

        type(integrator_config_t), intent(in) :: config

        if (config%integ_type == "rk45" .and. config%spatial_coords == "cyl") then
            integ = rk45_cyl_integrator_t(1d30, config%ro0, config%rtol)
        else if (config%integ_type == "rk45" .and. &
                 config%spatial_coords == "cyl_can") then
            integ = rk45_can_integrator_t(1d30, config%ro0, config%rtol)
        else
            print *, "create_integrator_vpar: Unknown type ", config%integ_type
            error stop
        end if

    end function create_integrator_vpar

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine integrator_init(si, field, f, config)

        type(symplectic_integrator_data_t), intent(inout) :: si
        class(field_can_t), intent(in) :: field
        type(field_can_data_t), intent(inout) :: f
        type(integrator_config_t), intent(in) :: config

        si%atol = 1d-15
        si%rtol = config%rtol

        si%ntau = config%nskip
        si%dt = config%dt

        si%z = config%zstart(1:4)

        call field%evaluate(f, si%z(1), si%z(2), si%z(3), 0)
        call get_val(f, si%z(4)) ! for pth
        si%pthold = f%pth
    end subroutine integrator_init


end module integrator
