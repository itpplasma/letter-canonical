program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use magfie, only: FieldType
    use magfie_factory, only: magfie_type_from_string
    use magfie_tok, only: input_file_tok => input_file
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components, init_splines_with_psi
    use field_can, only: field_can_t, field_can_cyl_t, field_can_new_t
    use integrator

    implicit none

    character(1024) :: input_file = "letter_canonical.in"

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer :: nt = 8000
    real(dp) :: ro0 = 1d0*20000d0  ! 1cm Larmor radius at 20000 Gauss
    real(dp) :: mu = 0d0 !1d-5

    class(FieldType), allocatable :: field_type
    class(field_can_t), allocatable :: field
    type(field_can_data_t) :: f
    class(symplectic_integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si

    real(dp) :: rmin, rmax, zmin, zmax
    real(dp) :: psi0 = 15000000d0
    real(dp) :: phi0 = 0d0
    real(dp) :: th0 = 27.5d0
    real(dp) :: vpar0 = 1.0d0
    real(dp) :: z0(4), starttime, endtime, dt = 1.0d0, rtol = 1d-13
    real(dp), allocatable :: out(:, :)

    integer :: kt, ierr

    ! Configuration in letter_canonical.in
    character(16) :: magfie_type = "tok"
    character(16) :: integrator_type = "euler1"
    character(1024) :: outname = "orbit.out"
    namelist /letter_canonical/ magfie_type, integrator_type, outname, input_file_tok, &
        ro0, mu, nt, dt, rtol, psi0, phi0, th0, vpar0

    call set_bounding_box
    call read_input_file

    field_type = magfie_type_from_string(trim(magfie_type))

    print *, "init_canonical ..."
    call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
        [rmax, twopi, zmax], field_type)

    print *, "init_transformation ..."
    call init_transformation

    print *, "init_canonical_field_components ..."
    call init_canonical_field_components

    print *, "init_splines_with_psi ..."
    call init_splines_with_psi

    call init_field_can
    call set_initial_conditions
    call init_integrator

    starttime = omp_get_wtime()
    call trace_orbit
    endtime = omp_get_wtime()

    print *, trim(outname), endtime-starttime

    call write_output

contains

    subroutine set_bounding_box
    ! Workaround, otherwise not initialized without perturbation field
        rmin = 75.d0
        rmax = 264.42281879194627d0
        zmin = -150.d0
        zmax = 147.38193979933115d0
    end subroutine set_bounding_box

    subroutine read_input_file
        integer :: iunit
        open(newunit=iunit, file=input_file, status='old')
        read(iunit, nml=letter_canonical)
        close(iunit)
    end subroutine read_input_file

    subroutine init_field_can
        field = field_can_new_t()
        call field_can_init(f, mu, ro0, vpar0)
        call field%evaluate(f, psi0, th0, phi0, 2)
    end subroutine init_field_can

    subroutine set_initial_conditions
        z0(1) = psi0   ! psi_c
        z0(2) = th0    ! Z
        z0(3) = phi0   ! varphi_c
        z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi

        allocate(out(5,nt))

        out=0d0

        out(1:4,1) = z0
        out(5,1) = f%H
    end subroutine set_initial_conditions

    subroutine init_integrator
        integ = create_integrator(trim(integrator_type), field)
        call integrator_init(si, field, f, z0, dt, 1, rtol)
    end subroutine init_integrator

    subroutine trace_orbit
        do kt = 2, nt
            ierr = 0
            call integ%timestep(si, f, ierr)
            if (.not. ierr==0) error stop
            out(1:4,kt) = si%z
            out(5,kt) = f%H
        end do
    end subroutine trace_orbit

    subroutine write_output
        open(unit=20, file=outname, action='write', recl=4096)
        do kt = 1, nt
            write(20,*) out(:,kt)
        end do
        close(20)
    end subroutine write_output

end program main
