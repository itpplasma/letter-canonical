module letter_canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_tok, only: tok_field_t
    use integrator, only: integrator_t
    use canonical, only: init_canonical, twopi
    implicit none

    integer :: error_code_ = 0

    integer, parameter :: n_r=100, n_phi=64, n_z=75

    class(tok_field_t), allocatable :: field
    class(integrator_t), allocatable :: integ

    real(dp) :: dtau_
    integer :: ntau_, nskip_

    abstract interface
        subroutine callback_type(z)
            import :: dp
            real(dp), intent(in) :: z(5)
        end subroutine callback_type
    end interface

contains

    subroutine init(field_file, integrator_type, dtau, ntau, nskip)
        use magfie_tok, only: input_file

        character(*), intent(in) :: field_file
        character(*), intent(in) :: integrator_type  ! rk45, expl_euler, expl_impl_euler
        ! Cylindrical coordinates
        ! Canonical coordinates with vpar as variable
        ! Canonical coordinates with pphi as variable
        real(dp), intent(in) :: dtau
        integer, intent(in) :: ntau, nskip

        real(dp) :: rmin, rmax, zmin, zmax

        dtau_ = dtau
        ntau_ = ntau
        nskip_ = nskip

        input_file = field_file
        call get_bounding_box(rmin, rmax, zmin, zmax)

        field = tok_field_t()

        print *, "init_magfie ..."
        call field%init_magfie()

        if (integrator_type == "rk45_cyl") then
            ! TODO: integ = rk45_cyl_integrator_t()
            print *, "using cylindrical coordinates"
            return
        end if

        print *, "using canonical cylindrical coordinates"
        print *, "init_canonical ..."

        call init_canonical( &
            n_r, n_phi, n_z, [rmin, 0.0d0, zmin], [rmax, twopi, zmax], field)

        call init_integrator_can(integrator_type)

        if (.not. allocated(integ)) then
            call throw_error("init: unknown integrator type")
            return
        end if

    end subroutine init


    subroutine init_integrator_can(integrator_type)
        character(*), intent(in) :: integrator_type

        if (integrator_type == "expl_euler") then
            ! TODO: integ = expl_euler_can_integrator_t()
            return
        else if (integrator_type == "rk45_can") then
            ! TODO: integ = rk45_can_integrator_t()
            return
        end if

        call init_integrator_pphi(integrator_type)
    end subroutine init_integrator_can


    subroutine init_integrator_pphi(integrator_type)
        character(*), intent(in) :: integrator_type

        if (integrator_type == "expl_impl_euler") then
            ! TODO: integ = expl_impl_euler_integrator_t()
        end if
    end subroutine init_integrator_pphi


    subroutine trace_orbit(z0, z_out, callback)
        real(dp), intent(in) :: z0(5)
        real(dp), intent(inout) :: z_out(:, :)
        procedure(callback_type), optional :: callback

        integer :: kt, ierr
        real(dp) :: z(5)

        if (.not. size(z_out, 2) == ntau_/nskip_) then
            call throw_error("trace_orbit: z_out has wrong size")
            return
        end if

        z_out(:, 1) = z0
        do kt = 2, ntau_
            call integ%timestep(z, dtau_, ierr)
            if (ierr /= 0) then
                call throw_error("trace_orbit: error in timestep", ierr)
                return
            end if
            call callback(z)
            if (mod(kt, nskip_) == 0) then
                z_out(:, kt/nskip_) = z
            end if
        end do
    end subroutine trace_orbit


    subroutine get_bounding_box(rmin, rmax, zmin, zmax)
        ! Workaround, otherwise not initialized without perturbation field
        real(dp), intent(out) :: rmin, rmax, zmin, zmax
        rmin = 75.d0
        rmax = 264.42281879194627d0
        zmin = -150.d0
        zmax = 147.38193979933115d0
    end subroutine get_bounding_box


    subroutine throw_error(msg, error_code)
        character(*), intent(in) :: msg
        integer, intent(in), optional :: error_code
        print *, msg
        error_code_ = error_code
    end subroutine throw_error

end module letter_canonical
