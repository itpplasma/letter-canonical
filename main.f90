program main
    use magfie_factory, only: magfie_type_from_string
    use magfie, only: FieldType
    use magfie_tok, only: TokFieldType
    use magfie_test, only: TestFieldType
    use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der2, destroy_splines_3d
    use canonical, only: init_canonical, init_transformation, twopi, &
        init_canonical_field_components, init_splines_with_psi, &
        spl_lam, spl_chi, spl_A1, spl_A2, spl_A3

    implicit none
    save

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer :: outfile_unit
    real(8) :: rmin, rmax, zmin, zmax
    !complex(8) :: pert

    class(FieldType), allocatable :: field_type

    ! Workaround, otherwise not initialized without perturbation field
    rmin = 75.d0
    rmax = 264.42281879194627d0
    zmin = -150.d0
    zmax = 147.38193979933115d0

    field_type = TokFieldType()
    !field_type = TestFieldType()
    ! pert = dcmplx(2.0d3, 2.0d3)
    ! call field_type%add_perturbation(3, -2, [pert, pert, pert])
    ! call field_type%add_perturbation(3, 2, [pert, pert, pert])
    ! call field_type%add_perturbation(5, -2, [pert, pert, pert])
    ! call field_type%add_perturbation(5, 2, [pert, pert, pert])

    print *, "init_canonical ..."
    call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
        [rmax, twopi, zmax], field_type)

    print *, "init_transformation ..."
    call init_transformation

    print *, "init_canonical_field_components ..."
    call init_canonical_field_components

    print *, "init_splines_with_psi ..."
    call init_splines_with_psi

    print *, "test_splines ..."
    call test_splines

    print *, "test_integration ..."
    call test_integration

    print *, "test_psi ..."
    call test_psi


contains

    subroutine test_psi
        use canonical, only: psi_of_x, psi_grid
        real(8), dimension(3) :: x
        real(8) :: psi
        integer :: i_r, i_phi, i_z

        i_r = 50
        i_phi = 32
        i_z = 37

        x = [rmin + (rmax - rmin) * (i_r - 1) / (n_r - 1), &
             twopi * (i_phi - 1) / (n_phi - 1), &
             zmin + (zmax - zmin) * (i_z - 1) / (n_z - 1)]

        psi = psi_of_x(i_r, i_phi, i_z)

        print *, "psi_of_x(i_r, i_phi, i_z) = ", psi_of_x(i_r, i_phi, i_z)
        print *, "psi_linspace(i_r) = ", psi_grid(i_r)

        deallocate(psi_of_x)
    end subroutine test_psi

    subroutine test_integration
        real(8), parameter :: tol = 1.0d-10
        real(8), parameter :: dt = 5.75d-3*twopi
        integer, parameter :: nt = 3000

        real(8) :: x0(3), x(3), xcyl(3), lam
        integer :: i_t, i_fs, n_flux

        n_flux = 15

        do i_fs = 1, n_flux
            x0 = [170.0d0, 0.0d0, 20.0d0]
            x0(3) = (i_fs*1.0d0/n_flux - 0.5d0)*130d0
            x = x0
            do i_t = 0, nt
                call odeint_allroutines(x, 3, i_t*dt, (i_t+1)*dt, tol, Bnoncan)
                write(100, *) x
            end do
        end do

        do i_fs = 1, n_flux
            x0 = [170.0d0, 0.0d0, 20.0d0]
            x0(3) = (i_fs*1.0d0/n_flux - 0.5d0)*130d0
            x = x0
            do i_t = 0, nt
                call odeint_allroutines(x, 3, i_t*dt, (i_t+1)*dt, tol, Bcan)
                call evaluate_splines_3d(spl_lam, x, lam)
                xcyl(1) = x(1)
                xcyl(2) = modulo(x(2) + lam, twopi)
                xcyl(3) = x(3)
                write(101, *) xcyl
                write(102, *) x
            end do
        end do
    end subroutine test_integration


    subroutine Bnoncan(t, x, dx)
        use canonical, only: magfie_type
        use magfie, only: compute_bfield

        real(8), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(3), intent(in) :: x
        real(8), dimension(3), intent(inout) :: dx
        real(8) :: BR, BZ, Bphi, Bphictr

        call magfie_type%compute_bfield(x(1), modulo(x(2), twopi), x(3), BR, Bphi, BZ)

        Bphictr = Bphi/x(1)  ! contravariant component
        dx(1) = BR/Bphictr
        dx(2) = 1d0
        dx(3) = BZ/Bphictr

    end subroutine Bnoncan


    subroutine Bcan(t, x, dx)
        use magfie_test, only: AMPL, AMPL2
        use canonical, only: evaluate_afield_can

        real(8), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(3), intent(in) :: x
        real(8), dimension(3), intent(inout) :: dx
        real(8) :: A1s, A2s, A3s, dA1s(3), dA2s(3), dA3s(3)
        real(8) :: B(3)

        call evaluate_afield_can(x, A1s, dA1s, A2s, dA2s, A3s, dA3s)

        B(1) = dA3s(2) - dA2s(3)
        B(2) = -dA3s(1)
        B(3) = dA2s(1)

        dx(1) = B(1)/B(2)
        dx(2) = 1d0
        dx(3) = B(3)/B(2)
    end subroutine Bcan


    subroutine test_splines

        integer, parameter :: n_r_test=49, n_phi_test=63, n_z_test=74

        real(8), dimension(:,:,:), allocatable :: lam_test, chi_test
        real(8), dimension(:,:,:), allocatable :: A1_test, A2_test, A3_test
        real(8), dimension(:,:,:,:), allocatable :: dlam_test, dchi_test
        real(8), dimension(3) :: x
        real(8) :: dummy(6)

        integer :: i_r, i_phi, i_z

        allocate(lam_test(n_r_test, n_phi_test, n_z_test))
        allocate(chi_test(n_r_test, n_phi_test, n_z_test))
        allocate(A1_test(n_r_test, n_phi_test, n_z_test))
        allocate(A2_test(n_r_test, n_phi_test, n_z_test))
        allocate(A3_test(n_r_test, n_phi_test, n_z_test))
        allocate(dlam_test(3, n_r_test, n_phi_test, n_z_test))
        allocate(dchi_test(3, n_r_test, n_phi_test, n_z_test))

        do i_z = 1, n_z_test
            do i_phi = 1, n_phi_test
                do i_r = 1, n_r_test
                    x = [rmin + (rmax - rmin) * (i_r - 1) / (n_r_test - 1), &
                         twopi * (i_phi - 1) / (n_phi_test - 1), &
                         zmin + (zmax - zmin) * (i_z - 1) / (n_z_test - 1)]
                    call evaluate_splines_3d_der2(spl_lam, x, &
                        lam_test(i_r, i_phi, i_z), dlam_test(:, i_r, i_phi, i_z), dummy)
                    call evaluate_splines_3d_der2(spl_chi, x, &
                        chi_test(i_r, i_phi, i_z), dchi_test(:, i_r, i_phi, i_z), dummy)
                    call evaluate_splines_3d(spl_A1, x, A1_test(i_r, i_phi, i_z))
                    call evaluate_splines_3d(spl_A2, x, A2_test(i_r, i_phi, i_z))
                    call evaluate_splines_3d(spl_A3, x, A3_test(i_r, i_phi, i_z))
                end do
            end do
        end do

        open(newunit=outfile_unit, file="lam_spl.out")
            write(outfile_unit, *) lam_test
        close(outfile_unit)

        open(newunit=outfile_unit, file="dlamdr_spl.out")
            write(outfile_unit, *) dlam_test(1,:,:,:)
        close(outfile_unit)

        open(newunit=outfile_unit, file="dlamdz_spl.out")
            write(outfile_unit, *) dlam_test(2,:,:,:)
        close(outfile_unit)

        open(newunit=outfile_unit, file="dlamdp_spl.out")
            write(outfile_unit, *) dlam_test(3,:,:,:)
        close(outfile_unit)


        open(newunit=outfile_unit, file="chi_spl.out")
            write(outfile_unit, *) chi_test
        close(outfile_unit)

        open(newunit=outfile_unit, file="dchidr_spl.out")
            write(outfile_unit, *) dchi_test(1,:,:,:)
        close(outfile_unit)

        open(newunit=outfile_unit, file="dchidz_spl.out")
            write(outfile_unit, *) dchi_test(2,:,:,:)
        close(outfile_unit)

        open(newunit=outfile_unit, file="dchidp_spl.out")
            write(outfile_unit, *) dchi_test(3,:,:,:)
        close(outfile_unit)


        open(newunit=outfile_unit, file="A1_spl.out")
            write(outfile_unit, *) A1_test
        close(outfile_unit)
        open(newunit=outfile_unit, file="A2_spl.out")
            write(outfile_unit, *) A2_test
        close(outfile_unit)
        open(newunit=outfile_unit, file="A3_spl.out")
            write(outfile_unit, *) A3_test
        close(outfile_unit)

        deallocate(lam_test, chi_test, A1_test, A2_test, A3_test)
    end subroutine test_splines

end program main
