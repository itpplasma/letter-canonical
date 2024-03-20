program test
    implicit none

    integer, parameter :: n_r=2, n_z=3, n_phi=4
    real(8) :: eps = 1d-10

    call setup
    call test_get_grid_point
    call test_generate_regular_grid
    call test_can_to_cyl
    call test_compute_Bmod
    call test_compute_hcan
    call test_compute_Acan

contains

    subroutine setup
        use canonical, only: init_canonical

        call print_test("Test Setup")
        call init_canonical(n_r, n_z, n_phi)
        call print_ok
    end subroutine setup


    subroutine test_get_grid_point
        use canonical, only: get_grid_point, pi
        use my_little_magfie, only : rmin, rmax, zmin, zmax

        real(8) :: x_expected(3), x_computed(3)

        call print_test("test_get_grid_point")

        x_computed = get_grid_point(1, 1, 1)
        x_expected = [rmin, zmin, 0.0d0]
        if (any(abs(x_computed - x_expected) > eps)) error stop

        x_computed = get_grid_point(n_r, n_z, n_phi)
        x_expected = [rmax, zmax, 2.0d0*pi]
        if (any(abs(x_computed - x_expected) > eps)) error stop

        call print_ok
    end subroutine test_get_grid_point


    subroutine test_generate_regular_grid
        use canonical, only: generate_regular_grid, pi
        use my_little_magfie, only : rmin, rmax, zmin, zmax

        real(8) :: grid_computed(3, n_r, n_z, n_phi)
        real(8) :: x_computed(3), x_expected(3)

        call print_test("test_generate_regular_grid")

        call generate_regular_grid(grid_computed)

        x_computed = grid_computed(:, 1, 1, 1)
        x_expected = [rmin, zmin, 0.0d0]
        if (any(abs(x_computed - x_expected) > eps)) error stop

        x_computed = grid_computed(:, n_r, n_z, n_phi)
        x_expected = [rmax, zmax, 2.0d0*pi]
        if (any(abs(x_computed - x_expected) > eps)) error stop

        call print_ok
    end subroutine test_generate_regular_grid

    subroutine test_compute_Bmod
        use canonical, only: compute_Bmod

        real(8), dimension(3, n_r, n_z, n_phi) :: B
        real(8), dimension(n_r, n_z, n_phi) :: Bmod_computed

        call print_test("test_compute_Bmod")

        B(:,:,:,:) = sqrt(1.0d0/3.0d0)

        call compute_Bmod(B, Bmod_computed)
        if (any(abs(Bmod_computed - 1.0d0) > eps)) error stop

        call print_ok
    end subroutine test_compute_Bmod


    subroutine test_compute_hcan
        use canonical, only: compute_hcan, compute_bmod, &
            generate_regular_grid, cyl_to_cov, spl_lam, spl_chi
        use interpolate, only: destroy_splines_3d

        real(8), dimension(3, n_r, n_z, n_phi) :: B, Bcov, x
        real(8), dimension(n_r, n_z, n_phi) :: Bmod
        real(8), dimension(2, n_r, n_z, n_phi) :: hcan_expected, hcan_computed

        call generate_regular_grid(x)

        B(:,:,:,:) = sqrt(1.0d0/3.0d0)
        Bcov = B
        call cyl_to_cov(x, Bcov)

        call compute_Bmod(B, Bmod)
        hcan_expected(1,:,:,:) = Bcov(2,:,:,:)/Bmod
        hcan_expected(2,:,:,:) = Bcov(3,:,:,:)/Bmod

        call construct_zero_spline(spl_lam)
        call construct_zero_spline(spl_chi)
        call compute_hcan(B, Bmod, hcan_computed)
        print *, hcan_computed(:,1,:,:)
        if (any(abs(hcan_computed - hcan_expected) > eps)) then
            call print_fail
            print *, "should match for identical transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)
        call destroy_splines_3d(spl_chi)

    end subroutine test_compute_hcan


    subroutine test_compute_Acan
        use canonical, only: compute_Acan, generate_regular_grid, cyl_to_cov, &
            spl_lam, spl_chi
        use interpolate, only: destroy_splines_3d

        real(8), dimension(3, n_r, n_z, n_phi) :: A, Acov, x
        real(8), dimension(2, n_r, n_z, n_phi) :: Acan_computed

        call print_test("test_compute_Acan")

        A(:,:,:,:) = 1.0d0
        Acov = A
        call generate_regular_grid(x)
        call cyl_to_cov(x, Acov)

        call construct_zero_spline(spl_lam)
        call construct_zero_spline(spl_chi)
        call compute_Acan(A, Acan_computed)
        if (any(abs(Acan_computed - Acov(2:3,:,:,:)) > eps)) then
            call print_fail
            print *, "should match for identical transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)
        call destroy_splines_3d(spl_chi)

        call construct_linear_spline(spl_lam)
        call construct_linear_spline(spl_chi)
        call compute_Acan(A, Acan_computed)
        if (all(abs(Acan_computed - Acov(2:3,:,:,:)) < eps)) then
            call print_fail
            print *, "must not match for linear transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)
        call destroy_splines_3d(spl_chi)

        call print_ok
    end subroutine test_compute_Acan


    subroutine test_can_to_cyl
        use canonical, only: can_to_cyl, spl_lam, generate_regular_grid
        use interpolate, only: destroy_splines_3d

        real(8), dimension(3, n_r, n_z, n_phi) :: xcan, xcyl_computed

        call generate_regular_grid(xcan)

        call print_test("test_can_to_cyl")
        call construct_zero_spline(spl_lam)
        call can_to_cyl(xcan, xcyl_computed)
        if (any(abs(xcyl_computed - xcan) > eps)) then
            call print_fail
            print *, "should match for identical transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)

        call construct_linear_spline(spl_lam)
        call can_to_cyl(xcan, xcyl_computed)
        if (all(abs(xcyl_computed - xcan) < eps)) then
            call print_fail
            print *, "must not match for linear transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)

        call print_ok
    end subroutine test_can_to_cyl


    subroutine construct_zero_spline(spl)
        use canonical, only: twopi
        use interpolate, only : SplineData3D, construct_splines_3d
        use my_little_magfie, only : rmin, rmax, zmin, zmax

        type(SplineData3D), intent(out) :: spl

        real(8) :: x_min(3), x_max(3)
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .False., .True.]

        real(8), dimension(n_r, n_z, n_phi) :: zeros

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]
        zeros(:, :, :) = 0.0d0

        call construct_splines_3d(x_min, x_max, zeros, order, periodic, spl)
    end subroutine construct_zero_spline


    subroutine construct_linear_spline(spl)
        use canonical, only: twopi
        use interpolate, only : SplineData3D, construct_splines_3d
        use my_little_magfie, only : rmin, rmax, zmin, zmax

        type(SplineData3D), intent(out) :: spl

        real(8) :: x_min(3), x_max(3)
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .False., .True.]

        real(8), dimension(n_r, n_z, n_phi) :: linear

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]

        call fill_linear(linear)

        call construct_splines_3d(x_min, x_max, linear, order, periodic, spl)
    end subroutine construct_linear_spline


    subroutine fill_linear(x)
        real(8), intent(inout) :: x(n_r, n_z, n_phi)

        integer :: i_r, i_z, i_phi

        do i_phi = 1, n_phi
            do i_z = 1, n_z
                do i_r = 1, n_r
                    x(i_r, i_z, i_phi) = i_r + i_z + i_phi + 1.0d0
                end do
            end do
        end do
    end subroutine fill_linear


    subroutine print_test(test_name)

        character(*) :: test_name
        print *, "==> ", test_name
    end subroutine print_test


    subroutine print_ok
        print *, "    .................................................... OK"
    end subroutine print_ok


    subroutine print_fail
        print *, "    .................................................... FAIL"
    end subroutine print_fail

end program test
