program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use magfie, only: FieldType
    use magfie_factory, only: magfie_type_from_string
    use magfie_tok, only: input_file_tok_ => input_file
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components, init_splines_with_psi
    use field_can, only: field_can_t, field_can_cyl_t, field_can_new_t
    use integrator

    implicit none

    character(1024) :: input_file = "letter_canonical.in"
    character(1024) :: input_file_tok = "field_divB0.inp"

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer :: nt = 80000
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
    real(dp) :: th0 = -40.0d0
    real(dp) :: vpar0 = 1.0d0
    real(dp) :: z0(4), starttime, endtime, dt = 1.0d0, rtol = 1d-13
    real(dp), allocatable :: out(:, :)

    integer, parameter :: nplagr=4, nder=0

    ! Configuration in letter_canonical.in
    character(16) :: magfie_type = "tok"
    character(16) :: integrator_type = "euler1"
    character(1024) :: outfile = "orbit.out"
    namelist /letter_canonical/ magfie_type, integrator_type, outfile, input_file_tok, &
        ro0, mu, nt, dt, rtol, psi0, phi0, th0, vpar0

    call set_bounding_box

    if (command_argument_count() > 0) then
        call get_command_argument(1, input_file)
    end if

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

    print *, trim(outfile), endtime-starttime

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

        input_file_tok_ = trim(input_file_tok)
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
        integer :: kt, ierr
        logical :: cut = .false.
        real(dp) :: zcut(4)
        do kt = 2, nt
            ierr = 0
            call integ%timestep(si, f, ierr)
            if (.not. ierr==0) error stop
            out(1:4, kt) = si%z
            out(5, kt) = f%H
            if (kt > nplagr/2) then
                cut = detect_cut(out, kt, zcut)
            end if
            if (cut) then
                write(666,*) zcut
            end if
        end do
    end subroutine trace_orbit

    function detect_cut(z, kt, zcut) result(cut)
        use plag_coeff_sub, only: plag_coeff

        real(dp), intent(in) :: z(:,:)
        integer, intent(in) :: kt
        real(dp), intent(out) :: zcut(4)
        logical :: cut


        real(dp) :: fper = twopi
        real(dp) :: coef(0:nder,nplagr)
        integer :: iper, kper_prev, kper, i

        iper = nplagr/2+1
        kper_prev=int(z(3, kt - nplagr/2)/fper)

        cut = .true.
        do i = kt - nplagr/2 + 1, kt
            kper = int(z(3, i)/fper)
            if (kper == kper_prev) then
                cut = .false.
                exit
            end if
        end do

        if (.not. cut) return

        print *, kper*fper, z(3, (kt-nplagr+1):kt)
        call plag_coeff(nplagr, nder, kper*fper, z(3, (kt-nplagr+1):kt), coef)
        zcut = matmul(z(:, (kt-nplagr+1):kt), coef(0,:))

    end function detect_cut

    subroutine write_output
        integer :: kt

        open(unit=20, file=outfile, action='write', recl=4096)
        do kt = 1, nt
            write(20,*) out(:,kt)
        end do
        close(20)
    end subroutine write_output

end program main
