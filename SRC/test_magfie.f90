program test_magfie
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_sub, only: odeint_allroutines
    use magfie_factory, only: magfie_type_from_string
    use magfie, only: field_t
    use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der2, destroy_splines_3d
    use canonical, only: init_canonical, init_transformation, twopi, &
        init_canonical_field_components, init_splines_with_psi, &
        spl_lam, spl_A3, spl_R_of_xc, spl_Aphi_of_xc

    implicit none
    save

    integer, parameter :: n_r=100, n_phi=16, n_z=75
    integer :: outfile_unit
    real(dp) :: rmin, rmax, zmin, zmax
    !complex(8) :: pert

    class(field_t), allocatable :: field_type

    character(16) :: magfie_type="tok"

    ! Workaround, otherwise not initialized without perturbation field
    rmin = 75.d0
    rmax = 264.42281879194627d0
    zmin = -150.d0
    zmax = 147.38193979933115d0

    field_type = magfie_type_from_string(trim(magfie_type))

    !field_type = tok_field_t()
    !field_type = test_field_t()
    ! pert = dcmplx(2.0d3, 2.0d3)
    ! call field_type%add_perturbation(3, -2, [pert, pert, pert])
    ! call field_type%add_perturbation(3, 2, [pert, pert, pert])
    ! call field_type%add_perturbation(5, -2, [pert, pert, pert])
    ! call field_type%add_perturbation(5, 2, [pert, pert, pert])

    print *, "init_magfie ..."
    call field_type%init_magfie()

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
        use canonical, only: psi_of_x, psi_grid, psi_inner, psi_outer, dpsi_dR_positive
        real(dp), dimension(3) :: x
        real(dp) :: psi, R
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

        print *, "psi_inner = ", psi_inner, "psi_outer = ", psi_outer
        print *, "h_psi = ", (psi_outer-psi_inner)/dble(n_R-1)
        print *, "dpsi_dR_positive = ", dpsi_dR_positive

        call write_out("psi_of_x.out", psi_of_x)

        call evaluate_splines_3d(spl_R_of_xc, [psi, x(2), x(3)], R)

        print *, "test_psi: R = ", x(1), "=", R

    end subroutine test_psi

    subroutine test_integration
        real(dp), parameter :: tol = 1.0d-10
        real(dp), parameter :: dt = 5.75d-3*twopi
        integer, parameter :: nt = 30000
        integer, parameter :: n_flux = 15

        call test_integration_noncan(tol, dt, nt, n_flux)
        call test_integration_can(tol, dt, nt, n_flux)
        call test_integration_can2(tol, dt, nt, n_flux)
    end subroutine test_integration


    subroutine test_integration_noncan(tol, dt, nt, n_flux)
        real(dp), intent(in) :: tol, dt
        integer, intent(in) :: nt, n_flux

        real(dp) :: x0(3), x(3), lam
        integer :: i_t, i_fs

        do i_fs = 1, n_flux
            x0 = [170.0d0, 0.0d0, 20.0d0]
            x0(3) = (i_fs*1.0d0/n_flux - 0.5d0)*130d0
            call evaluate_splines_3d(spl_lam, x0, lam)
            x0(2) = modulo(x0(2) + lam, twopi)
            x = x0
            do i_t = 0, nt
                call odeint_allroutines(x, 3, i_t*dt, (i_t+1)*dt, tol, Bnoncan)
                write(100, *) x
            end do
        end do
    end subroutine


    subroutine Bnoncan(t, x, dx)
        use canonical, only: magfie_type
        use magfie, only: compute_bfield

        real(dp), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(dp), dimension(3), intent(in) :: x
        real(dp), dimension(3), intent(inout) :: dx
        real(dp) :: BR, BZ, Bphi, Bphictr

        call magfie_type%compute_bfield(x(1), modulo(x(2), twopi), x(3), BR, Bphi, BZ)

        Bphictr = Bphi/x(1)  ! contravariant component
        dx(1) = BR/Bphictr
        dx(2) = 1d0
        dx(3) = BZ/Bphictr

    end subroutine Bnoncan


    subroutine test_integration_can(tol, dt, nt, n_flux)
        real(dp), intent(in) :: tol, dt
        integer, intent(in) :: nt, n_flux

        real(dp) :: x0(3), x(3), xcyl(3), lam
        integer :: i_t, i_fs

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
    end subroutine


    subroutine Bcan(t, x, dx)
        use canonical, only: evaluate_afield_can

        real(dp), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(dp), dimension(3), intent(in) :: x
        real(dp), dimension(3), intent(inout) :: dx
        real(dp) :: A1s, A2s, A3s, dA1s(3), dA2s(3), dA3s(3), d2A1s(6), &
                   d2A2s(6), d2A3s(6)
        real(dp) :: B(3)

        call evaluate_afield_can(x, A1s, dA1s, d2A1s, A2s, dA2s, d2A2s, &
            A3s, dA3s, d2A3s)

        B(1) = dA3s(2) - dA2s(3)
        B(2) = -dA3s(1)
        B(3) = dA2s(1)

        dx(1) = B(1)/B(2)
        dx(2) = 1d0
        dx(3) = B(3)/B(2)
    end subroutine Bcan


    subroutine test_integration_can2(tol, dt, nt, n_flux)
        real(dp), intent(in) :: tol, dt
        integer, intent(in) :: nt, n_flux

        real(dp) :: x0(3), xc(3), x(3), xcyl(3), lam, psi0
        integer :: i_t, i_fs

        do i_fs = 1, n_flux
            x0 = [170.0d0, 0.0d0, 20.0d0]
            x0(3) = (i_fs*1.0d0/n_flux - 0.5d0)*130d0
            xc(2:3) = x0(2:3)
            call evaluate_splines_3d(spl_A3, x0, psi0)
            xc(1) = psi0
            do i_t = 0, nt
                call odeint_allroutines(xc, 3, i_t*dt, (i_t+1)*dt, tol, Bcan2)
                x(2:3) = xc(2:3)
                call evaluate_splines_3d(spl_R_of_xc, xc, x(1))
                call evaluate_splines_3d(spl_lam, x, lam)
                xcyl(1) = x(1)
                xcyl(2) = modulo(x(2) + lam, twopi)
                xcyl(3) = x(3)
                write(103, *) xcyl
                write(104, *) x
                write(105, *) xc
            end do
        end do
    end subroutine


    subroutine Bcan2(t, xc, dxc)
        real(dp), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(dp), dimension(3), intent(in) :: xc
        real(dp), dimension(3), intent(inout) :: dxc
        real(dp) :: Aphi, dAphi(3), dummy(6)

        call evaluate_splines_3d_der2(spl_Aphi_of_xc, xc, Aphi, dAphi, dummy)

        dxc(1) = dAphi(3)
        dxc(2) = 1d0
        dxc(3) = -dAphi(1)
    end subroutine Bcan2


    subroutine test_splines
        use canonical

        integer, parameter :: n_r_test=49, n_phi_test=63, n_z_test=74

        real(dp), dimension(:,:,:), allocatable :: lam_test, chi_test
        real(dp), dimension(:,:,:), allocatable :: A2_test, A3_test
        real(dp), dimension(:,:,:,:), allocatable :: dlam_test, dchi_test
        real(dp), dimension(3) :: x
        real(dp) :: dummy(6)

        integer :: i_r, i_phi, i_z

        allocate(lam_test(n_r_test, n_phi_test, n_z_test))
        allocate(chi_test(n_r_test, n_phi_test, n_z_test))
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
                    call evaluate_splines_3d(spl_A2, x, A2_test(i_r, i_phi, i_z))
                    call evaluate_splines_3d(spl_A3, x, A3_test(i_r, i_phi, i_z))
                end do
            end do
        end do

        call write_out("lam_spl.out", lam_test)
        call write_out("dlamdr_spl.out", dlam_test(1,:,:,:))
        call write_out("dlamdp_spl.out", dlam_test(2,:,:,:))
        call write_out("dlamdz_spl.out", dlam_test(3,:,:,:))

        call write_out("chi_spl.out", chi_test)

        call write_out("dchidr_spl.out", dchi_test(1,:,:,:))
        call write_out("dchidp_spl.out", dchi_test(2,:,:,:))
        call write_out("dchidz_spl.out", dchi_test(3,:,:,:))
        call write_out("A2_spl.out", A2_test)
        call write_out("A3_spl.out", A3_test)

        deallocate(lam_test, chi_test, A2_test, A3_test)

        call write_out("R_of_xc.out", R_of_xc)
        call write_out("Aph_of_xc.out", Aph_of_xc)
        call write_out("hph_of_xc.out", hph_of_xc)
        call write_out("hth_of_xc.out", hth_of_xc)
        call write_out("Bmod_of_xc.out", Bmod_of_xc)

    end subroutine test_splines

    subroutine write_out(filename, data)
        character(len=*), intent(in) :: filename
        real(dp), dimension(:,:,:), intent(in) :: data

        open(newunit=outfile_unit, file=filename)
            write(outfile_unit, *) data
        close(outfile_unit)
    end subroutine write_out

end program test_magfie
