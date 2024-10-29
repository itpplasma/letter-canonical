module letter_canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_tok, only: tok_field_t, input_file_magfie_tok => input_file
    use integrator
    use callback, only: callback_pointer_t, cut_callback_t
    use canonical, only: init_canonical, init_transformation, &
        init_canonical_field_components, init_splines_with_psi, &
        can_psi_to_cyl, cyl_to_can_psi, twopi
    use field_can, only: field_can_t, field_can_albert_t, field_can_data_t, get_val
    use field_can_base, only: n_field_evaluations

    implicit none

    integer, parameter :: MAX_PATH_LENGTH = 1024

    integer :: error_code_ = 0

    character(MAX_PATH_LENGTH) :: input_file_tok
    character(MAX_PATH_LENGTH) :: output_prefix

    character(64) :: magfie_type = "tok"
    character(64) :: integrator_type = "rk45"
    character(64) :: spatial_coordinates = "cyl"
    character(64) :: velocity_coordinate = "vpar"

    integer :: n_r=100, n_phi=4, n_z=150

    real(dp) :: rmin = 75.d0, &
        rmax = 264.42281879194627d0, &
        zmin = -150.d0, &
        zmax = 147.38193979933115d0

    real(dp) :: ZPLANE_VALUE = -105.0d0

    class(tok_field_t), allocatable :: field
    class(field_can_t), allocatable :: field_can_
    class(integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si
    type(field_can_data_t) :: f

    real(dp) :: rmu=1d30, ro0=20000d0*2d0  ! 20000 Gauss, 1cm Larmor radius

    real(dp) :: dtau=1d0, rtol=1d-13
    integer :: ntau=1000, nskip=1

    real(dp) :: R0=154.0d0, phi0=-6.283d0, Z0=0.0d0, vpar0=0d0

    real(dp) :: mu  ! magnetic moment

    real(dp), allocatable :: z_out(:, :)

    class(callback_pointer_t), allocatable :: callbacks(:)
    class(cut_callback_t), allocatable :: phi0_plus_callback, phi0_minus_callback, &
        vpar0_plus_callback, vpar0_minus_callback, zplane_callback

    integer :: phi0plus_unit, phi0minus_unit, vpar0plus_unit, vpar0minus_unit, &
        zplane_unit

    namelist /config/ magfie_type, integrator_type, input_file_tok, &
        output_prefix, spatial_coordinates, velocity_coordinate, rtol, &
        rmin, rmax, zmin, zmax, rmu, ro0, dtau, ntau, nskip, n_r, n_phi, n_z, &
        R0, phi0, Z0, vpar0

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

        input_file_magfie_tok = input_file_tok

        print *, "init: main output is ", trim(output_prefix) // ".out"

        field = tok_field_t()

        print *, "init_magfie ..."
        call field%init_magfie()

        if (spatial_coordinates == "cyl") then
            print *, "init: using cylindrical coordinates"
        else if (spatial_coordinates == "albert") then
            print *, "init: using canonical cylindrical coordinates with A_Z as radius"
            call init_canonical(n_r, n_phi, n_z, &
                [rmin, 0.0d0, zmin], [rmax, twopi, zmax], field)
            call init_transformation
            call init_canonical_field_components
            call init_splines_with_psi
        else
            call throw_error("init: " // integ_error_message())
            return
        end if

        call init_integrator
        call init_callbacks

        allocate(z_out(5, ntau/nskip))

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

        if (i == 0) then
            prefix = path
        else
            prefix = path(1:i-1)
        end if
    end subroutine remove_extension


    subroutine init_integrator
        type(integrator_config_t) :: integ_config
        real(dp) :: z_internal(5)

        if (velocity_coordinate == "pphi") then
            field_can_ = create_field_can(spatial_coordinates)
        end if

        print *, "init_integrator ...."

        mu = compute_mu([R0, phi0, Z0, 1d0, vpar0])

        call to_internal_coordinates([R0, phi0, Z0, 1d0, vpar0], z_internal)

        integ_config = integrator_config_t(integrator_type, spatial_coordinates, &
            velocity_coordinate, z_internal, dtau, ro0, rtol, 1)

        if (integ_config%momentum_coord == "pphi") then
            integ = create_integrator(integ_config, si, field_can_, f)
        else
            integ = create_integrator(integ_config)
        end if
    end subroutine init_integrator


    function compute_mu(z)
        real(dp) :: compute_mu
        real(dp), intent(in) :: z(5)
        real(dp) :: BR, Bphi, BZ, Bmod

        call field%compute_bfield(z(1), z(2), z(3), BR, Bphi, BZ)

        Bmod = dsqrt(BR**2 + Bphi**2 + BZ**2)

        compute_mu = 0.5d0*z(4)**2*(1.d0-z(5)**2)/Bmod

    end function compute_mu


    subroutine init_callbacks
        integer, parameter :: NUM_CALLBACKS = 5
        allocate(callbacks(NUM_CALLBACKS))

        allocate(phi0_plus_callback)
        phi0_plus_callback%distance => phi0_plus_distance
        phi0_plus_callback%event => phi0_plus_write
        allocate(callbacks(1)%item, source=phi0_plus_callback)
        open(newunit=phi0plus_unit, file=trim(output_prefix) // "_phi0plus.out")

        allocate(phi0_minus_callback)
        phi0_minus_callback%distance => phi0_minus_distance
        phi0_minus_callback%event => phi0_minus_write
        allocate(callbacks(2)%item, source=phi0_minus_callback)
        open(newunit=phi0minus_unit, file=trim(output_prefix) // "_phi0minus.out")

        allocate(vpar0_plus_callback)
        vpar0_plus_callback%distance => vpar0_plus_distance
        vpar0_plus_callback%event => vpar0_plus_write
        allocate(callbacks(3)%item, source=vpar0_plus_callback)
        open(newunit=vpar0plus_unit, file=trim(output_prefix) // "_vpar0plus.out")

        allocate(vpar0_minus_callback)
        vpar0_minus_callback%distance => vpar0_minus_distance
        vpar0_minus_callback%event => vpar0_minus_write
        allocate(callbacks(4)%item, source=vpar0_minus_callback)
        open(newunit=vpar0minus_unit, file=trim(output_prefix) // "_vpar0minus.out")

        allocate(zplane_callback)
        zplane_callback%distance => zplane_distance
        zplane_callback%event => zplane_write
        allocate(callbacks(5)%item, source=zplane_callback)
        open(newunit=zplane_unit, file=trim(output_prefix) // "_zplane.out")
    end subroutine init_callbacks


    subroutine trace_orbit
        integer :: i, kt, ierr
        real(dp) :: zstart(5), z(5), zcyl(5)

        zstart = [R0, phi0, Z0, 1d0, vpar0]

        call to_internal_coordinates(zstart, z)

        z_out(:, 1) = z
        do kt = 2, ntau
            call integ%timestep(z, dtau, ierr)
            if (ierr /= 0) then
                call throw_error("trace_orbit: error in timestep", ierr)
                return
            end if
            call from_internal_coordinates(z, zcyl)
            do i = 1, size(callbacks)
                call callbacks(i)%execute(kt*dtau, zcyl)
            end do
            if (mod(kt, nskip) == 0) then
                z_out(:, kt/nskip) = z
            end if
        end do

        print *, "timesteps: ", ntau
        print *, "field evaluations: ", integ%get_field_evaluations()
    end subroutine trace_orbit


    subroutine to_internal_coordinates(z, z_internal)
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: z_internal(5)

        real(dp) :: x(3)

        if (spatial_coordinates == "cyl") then
            z_internal(1:3) = z(1:3)
        else if (spatial_coordinates == "albert") then
            call cyl_to_can_psi(z(1:3), x)
            z_internal(1:3) = x([1,3,2])  ! swap phi and theta order
        else
            call throw_error("to_internal_coordinates: " // integ_error_message())
            return
        end if

        if (velocity_coordinate == "vpar") then
            z_internal(4:5) = z(4:5)
        else if (velocity_coordinate == "pphi") then
            call field_can_%evaluate(f, z_internal(1), z_internal(2), z_internal(3), 0)
            n_field_evaluations = n_field_evaluations - 1  ! don't count for benchmark

            ! normalization of thermal velocity different by factor sqrt(2) - see docs
            f%mu = .5d0*z(4)**2*(1.d0-z(5)**2)/f%Bmod*2d0
            f%ro0 = ro0/dsqrt(2d0)
            f%vpar = z(4)*z(5)*dsqrt(2d0)

            ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
            z_internal(4) = f%vpar*f%hph + f%Aph/f%ro0

            ! we use our free last variable for the Hamiltonian H
            z_internal(5) = f%vpar**2/2d0 + f%mu*f%Bmod
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
        else if (spatial_coordinates == "albert") then
            call can_psi_to_cyl(z_internal([1,3,2]), z(1:3))  ! swap phi and theta order
        else
            call throw_error("from_internal_coordinates: " // integ_error_message())
            return
        end if

        if (velocity_coordinate == "vpar") then
            z(4:5) = z_internal(4:5)
        else if (velocity_coordinate == "pphi") then
            call field_can_%evaluate(f, z_internal(1), z_internal(2), z_internal(3), 0)
            n_field_evaluations = n_field_evaluations - 1  ! don't count for benchmark

            call get_val(f, z_internal(4))
            z(4) = z_internal(5)
            z(5) = f%vpar/(z_internal(5)*dsqrt(2d0))
        else
            call throw_error("from_internal_coordinates: " // integ_error_message())
            return
        end if
    end subroutine from_internal_coordinates


    subroutine write_output
        real(dp) :: z(5)
        integer :: kt, iunit, iunit_pphi_H
        character(MAX_PATH_LENGTH) :: outname, outname_pphi_H

        real(dp) :: pphi, H

        outname = trim(output_prefix) // ".out"
        outname_pphi_H = trim(output_prefix) // "_pphi_H.out"
        open(newunit=iunit, file=outname)
        open(newunit=iunit_pphi_H, file=outname_pphi_H)
        do kt = 1, size(z_out, 2)
            call from_internal_coordinates(z_out(:, kt), z)
            write(iunit, *) (kt-1)*dtau*nskip, z
            call get_pphi_H(z, pphi, H)
            write(iunit_pphi_H, *) (kt-1)*dtau*nskip, pphi, H
        end do
        close(iunit_pphi_H)
        close(iunit)

        outname = trim(output_prefix) // "_internal.out"
        open(newunit=iunit, file=outname)
        do kt = 1, size(z_out, 2)
            write(iunit, *) (kt-1)*dtau*nskip, z_out(:, kt)
        end do
        close(iunit)
    end subroutine write_output


    subroutine get_pphi_H(z, pphi, H)
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: pphi, H

        real(dp) :: AR, Aphi, AZ, BR, Bphi, BZ, Bmod

        call field%compute_abfield(z(1), z(2), z(3), AR, Aphi, AZ, BR, Bphi, BZ)

        Bmod = dsqrt(BR**2 + Bphi**2 + BZ**2)

        pphi = (z(4)*z(5)*Bphi/Bmod + Aphi/ro0)*z(1)
        H = 2d0*(z(4)**2*0.5d0*z(5)**2 + mu*Bmod)  ! energy normalization factor 2
    end subroutine get_pphi_H


    function phi0_plus_distance(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = modulo(z(2), twopi) - 0.5d0*twopi
    end function phi0_plus_distance

    subroutine phi0_plus_write(t, z)
        real(dp), intent(in) :: t, z(:)
        write(phi0plus_unit, *) t, z
    end subroutine phi0_plus_write


    function phi0_minus_distance(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = -modulo(z(2), twopi) + 0.5d0*twopi
    end function phi0_minus_distance

    subroutine phi0_minus_write(t, z)
        real(dp), intent(in) :: t, z(:)
        write(phi0minus_unit, *) t, z
    end subroutine phi0_minus_write


    function vpar0_plus_distance(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = z(5)
    end function vpar0_plus_distance

    subroutine vpar0_plus_write(t, z)
        real(dp), intent(in) :: t, z(:)
        write(vpar0plus_unit, *) t, z
    end subroutine vpar0_plus_write


    function vpar0_minus_distance(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = -z(5)
    end function vpar0_minus_distance

    subroutine vpar0_minus_write(t, z)
        real(dp), intent(in) :: t, z(:)
        write(vpar0minus_unit, *) t, z
    end subroutine vpar0_minus_write


    function zplane_distance(t, z) result(distance)
        real(dp) :: distance
        real(dp), intent(in) :: t, z(:)
        distance = -(z(3) - ZPLANE_VALUE)
    end function zplane_distance



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
