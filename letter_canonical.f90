module letter_canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_tok, only: tok_field_t
    use integrator, only: integrator_t, rk45_cyl_integrator_t
    use canonical, only: init_canonical, twopi
    implicit none

    integer, parameter :: MAX_PATH_LENGTH = 1024

    integer :: error_code_ = 0

    character(MAX_PATH_LENGTH) :: input_file_tok
    character(MAX_PATH_LENGTH) :: output_prefix

    character(64) :: magfie_type = "tok"
    character(64) :: integrator_type = "rk45"
    character(64) :: spatial_coordinates = "cyl"
    character(64) :: velocity_coordinate = "vpar"

    integer :: n_r=100, n_phi=64, n_z=75

    real(dp) :: rmin, rmax, zmin, zmax

    class(tok_field_t), allocatable :: field
    class(integrator_t), allocatable :: integ

    real(dp) :: dtau
    integer :: ntau, nskip

    namelist /config/ magfie_type, integrator_type, input_file_tok, &
        output_prefix, spatial_coordinates, velocity_coordinate, rmin, rmax, zmin, zmax

    abstract interface
        subroutine callback_type(z)
            import :: dp
            real(dp), intent(in) :: z(5)
        end subroutine callback_type
    end interface

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

        print *, "init: main output is ", trim(output_prefix) // ".out"

        field = tok_field_t()

        print *, "init_magfie ..."
        call field%init_magfie()

        if (spatial_coordinates == "cyl") then
            print *, "init: using cylindrical coordinates"
            call init_integrator_cyl
        else if (spatial_coordinates == "can") then
            print *, "init: using canonical coordinates"
            print *, "init_canonical ..."
            call init_canonical( n_r, n_phi, n_z, &
                [rmin, 0.0d0, zmin], [rmax, twopi, zmax], field)
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

        open(newunit=iunit, file=input_file)
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
            integ = rk45_cyl_integrator_t()
        else
            call throw_error("init_integrator_cyl: " // trim(integ_error_message()))
            return
        end if
    end subroutine init_integrator_cyl


    subroutine init_integrator_can

        if (velocity_coordinate == "vpar" .and. integrator_type == "rk45") then
                ! TODO: integ = rk45_can_vpar_integrator_t()
        else if (velocity_coordinate == "pphi" .and. &
                integrator_type == "expl_impl_euler") then
                ! TODO: integ = expl_impl_euler_integrator_t()
        else
            call throw_error("init_integrator_can: " // trim(integ_error_message()))
            return
        end if
    end subroutine init_integrator_can


    subroutine trace_orbit(z0, z_out, callback)
        real(dp), intent(in) :: z0(5)
        real(dp), intent(inout) :: z_out(:, :)
        procedure(callback_type), optional :: callback

        integer :: kt, ierr
        real(dp) :: z(5)

        if (.not. size(z_out, 2) == ntau/nskip) then
            call throw_error("trace_orbit: z_out has wrong size")
            return
        end if

        z_out(:, 1) = z0
        do kt = 2, ntau
            call integ%timestep(z, dtau, ierr)
            if (ierr /= 0) then
                call throw_error("trace_orbit: error in timestep", ierr)
                return
            end if
            call callback(z)
            if (mod(kt, nskip) == 0) then
                z_out(:, kt/nskip) = z
            end if
        end do
    end subroutine trace_orbit


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
