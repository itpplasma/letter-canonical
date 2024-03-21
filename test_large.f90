program test_large
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    integer, parameter :: n_r=100, n_z=75, n_phi=64

    call setup
    call test_physical_field

contains

    subroutine setup
        use canonical, only: init_canonical, twopi

        real(8) :: rmin, rmax, zmin, zmax

        call print_test("Test Setup")
        rmin = 75.d0
        rmax = 264.42281879194627d0
        zmin = -150.d0
        zmax = 147.38193979933115d0
        call init_canonical(n_r, n_z, n_phi, [rmin, zmin, 0.0d0], &
            [rmax, zmax, twopi])
        call print_ok
    end subroutine setup


    subroutine test_physical_field
        use canonical, only: get_physical_field, generate_regular_grid, &
            x_min, x_max, order, periodic
        use interpolate, only: SplineData3D, construct_splines_3d, &
            evaluate_splines_3d, evaluate_splines_3d_der2

        real(8), parameter :: tol=5.0d-3
        real(8), dimension(:,:,:,:), allocatable :: xcyl, B, Acyl
        type(SplineData3D) :: spl_AR, spl_AZ, spl_Aphicov, &
            spl_BR, spl_BZ, spl_Bphi

        real(8) :: BR, BZ, Bphi, AR, AZ, Aphicov, dAR(3), dAZ(3), dAphicov(3), &
            B_pitch(2), dummy(6)
        real(8) :: BR_expected, BZ_expected, Bphi_expected, B_pitch_expected(2)
        real(8) :: x(3), gii(3), sqrtg

        allocate(xcyl(3,n_r,n_z,n_phi))
        allocate(B(3,n_r,n_z,n_phi), Acyl(3,n_r,n_z,n_phi))

        call print_test("test_physical_field")

        call generate_regular_grid(xcyl)
        call get_physical_field(xcyl, B, Acyl)

        call construct_splines_3d(x_min, x_max, &
            B(1,:,:,:), order, periodic, spl_BR)
        call construct_splines_3d(x_min, x_max, &
            B(2,:,:,:), order, periodic, spl_BZ)
        call construct_splines_3d(x_min, x_max, &
            B(3,:,:,:), order, periodic, spl_Bphi)
        deallocate(B)

        call construct_splines_3d(x_min, x_max, &
            Acyl(1,:,:,:), order, periodic, spl_AR)
        call construct_splines_3d(x_min, x_max, &
            Acyl(2,:,:,:), order, periodic, spl_AZ)
        call construct_splines_3d(x_min, x_max, &
            Acyl(3,:,:,:)*xcyl(1,:,:,:), order, periodic, spl_Aphicov)
        deallocate(Acyl)

        x = 0.5d0*(x_min + x_max)
        call evaluate_splines_3d(spl_BR, x, BR)
        call evaluate_splines_3d(spl_BZ, x, BZ)
        call evaluate_splines_3d(spl_Bphi, x, Bphi)
        call evaluate_splines_3d_der2(spl_AR, x, AR, dAR, dummy)
        call evaluate_splines_3d_der2(spl_AZ, x, AZ, dAZ, dummy)
        call evaluate_splines_3d_der2(spl_Aphicov, x, Aphicov, dAphicov, dummy)

        gii = [1.d0, 1.d0, x(1)**2]
        sqrtg = sqrt(gii(1)*gii(2)*gii(3))

        BR_expected = (dAphicov(2) - dAZ(3))/sqrtg*sqrt(gii(1))
        BZ_expected = (dAR(3) - dAphicov(1))/sqrtg*sqrt(gii(2))
        Bphi_expected = (dAZ(1) - dAR(2))/sqrtg*sqrt(gii(3))

        B_pitch_expected(1) = BR_expected/Bphi_expected
        B_pitch_expected(2) = BZ_expected/Bphi_expected
        B_pitch = [BR/Bphi, BZ/Bphi]

        if (any(abs((B_pitch - B_pitch_expected)/maxval(B_pitch)) > tol)) then
            call print_fail
            print *, "B_pitch = ", B_pitch
            print *, "B_pitch_expected = ", B_pitch_expected
        else
            call print_ok
        end if


        call print_ok
    end subroutine test_physical_field


end program test_large
