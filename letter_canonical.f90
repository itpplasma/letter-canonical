module letter_canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_tok
    use integrator
    use callback
    use canonical, only: init_canonical, init_transformation, &
        init_canonical_field_components, init_splines_with_psi, &
        can_psi_to_cyl, cyl_to_can_psi, twopi
    implicit none

    integer, parameter :: MAX_PATH_LENGTH = 1024

    integer :: error_code_ = 0

    character(MAX_PATH_LENGTH) :: input_file_tok
    character(MAX_PATH_LENGTH) :: output_prefix

    character(64) :: magfie_type = "tok"
    character(64) :: integrator_type = "rk45"
    character(64) :: spatial_coordinates = "cyl"
    character(64) :: velocity_coordinate = "vpar"

    integer :: n_r=100, n_phi=8, n_z=150

    real(dp) :: rmin = 75.d0, &
        rmax = 264.42281879194627d0, &
        zmin = -150.d0, &
        zmax = 147.38193979933115d0

    class(tok_field_t), allocatable :: field
    class(integrator_t), allocatable :: integ

    real(dp) :: rmu=1d10, ro0=20000d0*2d0  ! 20000 Gauss, 1cm Larmor radius

    real(dp) :: dtau=1d0
    integer :: ntau=1000, nskip=1

    namelist /config/ magfie_type, integrator_type, input_file_tok, &
        output_prefix, spatial_coordinates, velocity_coordinate, &
        rmin, rmax, zmin, zmax, rmu, ro0, dtau, ntau, nskip, n_r, n_phi, n_z

contains

    subroutine init(input_file)
        character(MAX_PATH_LENGTH), intent(in), optional :: input_file

        if (present(input_file)) then
            call read_input_file(input_file)
            call remove_extension(input_file, output_prefix)
        else
            call read_input_file("letter_canonical.in")
            output_prefix = "letter_canonical"
        end if

        if (magfie_type /= "tok") then
            call throw_error("init: magfie_type must be 'tok'")
            return
        end if

        print *, "init: main output is ", trim(output_prefix) // ".out"

        field = tok_field_t()

        print *, "init_magfie ..."
        call field%init_magfie()

        if (spatial_coordinates == "cyl") then
            print *, "init: using cylindrical coordinates"
            call init_integrator_cyl
            return
        else if (spatial_coordinates == "cyl_can") then
            print *, "init: using canonical cylindrical coordinates"
            print *, "init_canonical ..."
            call init_canonical(n_r, n_phi, n_z, &
                [rmin, 0.0d0, zmin], [rmax, twopi, zmax], field)
            call init_transformation
            call init_canonical_field_components
            call init_splines_with_psi
        end if

        call init_integrator_can

        if (.not. allocated(integ)) then
            call throw_error("init: " // trim(integ_error_message()))
            return
        end if

    end subroutine init


    subroutine read_input_file(input_file)
        character(*), intent(in) :: input_file
        integer :: iunit

        open(newunit=iunit, file=input_file, status="old")
        read(iunit, nml=config)
        close(iunit)
    end subroutine read_input_file


    subroutine remove_extension(path, prefix)
        character(*), intent(in) :: path
        character(*), intent(inout) :: prefix
        integer :: i

        i = index(path, ".", back=.true.)

        print *, "remove_extension: i = ", i

        if (i == 0) then
            prefix = path
        else
            prefix = path(1:i-1)
        end if
    end subroutine remove_extension


    subroutine init_integrator_cyl
        if (integrator_type == "rk45" .and. velocity_coordinate == "vpar") then
            integ = rk45_cyl_integrator_t(rmu, ro0, 1d-8)
        else
            call throw_error("init_integrator_cyl: " // trim(integ_error_message()))
            return
        end if
    end subroutine init_integrator_cyl


    subroutine init_integrator_can
        if (velocity_coordinate == "vpar" .and. integrator_type == "rk45") then
            integ = rk45_can_integrator_t(rmu, ro0, 1d-8)
        else if (&
            velocity_coordinate == "pphi" .and. integrator_type=="expl_impl_euler") then
            integ = expl_impl_euler_integrator_t()
        else
            call throw_error("init_integrator_can: " // trim(integ_error_message()))
            return
        end if
    end subroutine init_integrator_can


    subroutine trace_orbit(z0, z_out, callbacks)
        real(dp), intent(in) :: z0(5)
        real(dp), allocatable, intent(out) :: z_out(:, :)
        class(callback_pointer_t), intent(inout), optional :: callbacks(:)

        integer :: i, kt, ierr
        real(dp) :: z(5), zcyl(5)

        allocate(z_out(5, ntau/nskip))

        call to_internal_coordinates(z0, z)

        z_out(:, 1) = z
        do kt = 2, ntau
            call integ%timestep(z, dtau, ierr)
            if (ierr /= 0) then
                call throw_error("trace_orbit: error in timestep", ierr)
                return
            end if
            if (present(callbacks)) then
                do i = 1, size(callbacks)
                    call to_internal_coordinates(z, zcyl)
                    call callbacks(i)%execute(kt*dtau, zcyl)
                end do
            end if
            if (mod(kt, nskip) == 0) then
                z_out(:, kt/nskip) = z
            end if
        end do
    end subroutine trace_orbit


    function get_pphi(z) result(pphi)
        real(dp), intent(in) :: z(5)
        real(dp) :: pphi

        pphi = z(4)
    end function get_pphi


    subroutine to_internal_coordinates(z, z_internal)
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: z_internal(5)

        real(dp) :: x(3)

        if (spatial_coordinates == "cyl") then
            z_internal(1:3) = z(1:3)
        else if (spatial_coordinates == "cyl_can") then
            call cyl_to_can_psi(z(1:3), x)
            z_internal(1:3) = x([1,3,2])  ! swap phi and theta order
        else
            call throw_error("to_internal_coordinates: " // integ_error_message())
            return
        end if

        if (velocity_coordinate == "vpar") then
            z_internal(4:5) = z(4:5)
        else
            call throw_error("to_internal_coordinates: " // integ_error_message())
            return
        end if
    end subroutine to_internal_coordinates


    subroutine from_internal_coordinates(z_internal, z)
        real(dp), intent(in) :: z_internal(5)
        real(dp), intent(out) :: z(5)

        if (spatial_coordinates == "cyl") then
            z(1:3) = z_internal(1:3)
        else if (spatial_coordinates == "cyl_can") then
            call can_psi_to_cyl(z_internal([1,3,2]), z(1:3))  ! swap phi and theta order
        else
            call throw_error("from_internal_coordinates: " // integ_error_message())
            return
        end if

        if (velocity_coordinate == "vpar") then
            z(4:5) = z_internal(4:5)
        else
            call throw_error("from_internal_coordinates: " // integ_error_message())
            return
        end if
    end subroutine from_internal_coordinates


    subroutine write_output(z_out)
        real(dp), intent(in) :: z_out(:, :)
        real(dp) :: z(5)
        integer :: kt, iunit
        character(MAX_PATH_LENGTH) :: outname

        outname = trim(output_prefix) // ".out"
        open(newunit=iunit, file=outname)
        do kt = 1, size(z_out, 2)
            call from_internal_coordinates(z_out(:, kt), z)
            write(iunit, *) z
        end do
        close(iunit)
    end subroutine write_output


    function integ_error_message() result(msg)
        character(1024) :: msg
        msg = "combination [" // &
            trim(magfie_type) // ", " // &
            trim(spatial_coordinates) // ", " // &
            trim(integrator_type) // ", " // &
            trim(velocity_coordinate) // "] not implemented"
    end function integ_error_message


    subroutine throw_error(msg, error_code)
        character(*), intent(in) :: msg
        integer, intent(in), optional :: error_code
        if (present(error_code)) then
            error_code_ = error_code
        else
            error_code_ = -1
        end if
        print *, "ERROR ", msg
    end subroutine throw_error


    subroutine stop_on_error
        if (error_code_ /= 0) then
            error stop error_code_
        end if
    end subroutine stop_on_error

end module letter_canonical
