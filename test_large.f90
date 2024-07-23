program test_large
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    integer, parameter :: n_r=100, n_phi=64, n_z=75

    call setup
    call test_physical_field

contains

    subroutine setup
        use magfie_test, only: test_field_t
        use canonical, only: init_canonical, twopi

        real(dp) :: rmin, rmax, zmin, zmax

        call print_test("Test Setup")
        rmin = 75.d0
        rmax = 264.42281879194627d0
        zmin = -150.d0
        zmax = 147.38193979933115d0
        call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
            [rmax, twopi, zmax], test_field_t())
        call print_ok
    end subroutine setup


    subroutine test_physical_field
        use canonical, only: get_physical_field, generate_regular_grid, &
            x_min, x_max, order, periodic
        use interpolate, only: SplineData3D, construct_splines_3d, &
            evaluate_splines_3d, evaluate_splines_3d_der2

        real(dp), parameter :: tol=5.0d-3
        real(dp), dimension(:,:,:,:), allocatable :: xcyl, B, Acyl
        type(SplineData3D) :: spl_AR, spl_AZ, spl_Aphicov, &
            spl_BR, spl_BZ, spl_Bphi

        real(dp) :: BR, BZ, Bphi, AR, AZ, Aphicov, dAR(3), dAZ(3), dAphicov(3), &
            B_pitch(2), dummy(6)
        real(dp) :: BR_expected, BZ_expected, Bphi_expected, B_pitch_expected(2)
        real(dp) :: x(3), gii(3), sqrtg

        allocate(xcyl(3, n_r, n_phi, n_z))
        allocate(B(3, n_r, n_phi, n_z), Acyl(3, n_r, n_phi, n_z))

        call print_test("test_physical_field")

        call generate_regular_grid(xcyl)
        call get_physical_field(xcyl, B, Acyl)

        call construct_splines_3d(x_min, x_max, B(1,:,:,:), order, periodic, spl_BR)
        call construct_splines_3d(x_min, x_max, B(2,:,:,:), order, periodic, spl_Bphi)
        call construct_splines_3d(x_min, x_max, B(3,:,:,:), order, periodic, spl_BZ)
        deallocate(B)

        call construct_splines_3d(x_min, x_max, Acyl(1,:,:,:), order, periodic, spl_AR)
        call construct_splines_3d(x_min, x_max, Acyl(2,:,:,:)*xcyl(1,:,:,:), order, &
            periodic, spl_Aphicov)
        call construct_splines_3d(x_min, x_max, Acyl(3,:,:,:), order, periodic, spl_AZ)
        deallocate(Acyl)

        x = 0.5d0*(x_min + x_max)
        x(2) = 0.2d0
        call evaluate_splines_3d(spl_BR, x, BR)
        call evaluate_splines_3d(spl_BZ, x, BZ)
        call evaluate_splines_3d(spl_Bphi, x, Bphi)
        call evaluate_splines_3d_der2(spl_AR, x, AR, dAR, dummy)
        call evaluate_splines_3d_der2(spl_AZ, x, AZ, dAZ, dummy)
        call evaluate_splines_3d_der2(spl_Aphicov, x, Aphicov, dAphicov, dummy)

        gii = [1.d0, x(1)**2, 1.d0]
        sqrtg = sqrt(gii(1)*gii(2)*gii(3))

        BR_expected = (dAZ(2) - dAphicov(3))/sqrtg*sqrt(gii(1))
        Bphi_expected = (dAR(2) - dAZ(1))/sqrtg*sqrt(gii(2))
        BZ_expected = (dAphicov(1) - dAR(3))/sqrtg*sqrt(gii(3))

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
