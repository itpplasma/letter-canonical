program test
    implicit none

    integer, parameter :: n_r=2, n_z=3, n_phi=4
    real(8) :: eps = 1d-10

    call setup
    call test_get_grid_point
    call test_generate_regular_grid
    call test_compute_A2can

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


    subroutine test_compute_A2can
        use canonical, only: compute_A2can, spl_lam, spl_chi
        use interpolate, only: destroy_splines_3d

        real(8), dimension(3, n_r, n_z, n_phi) :: A
        real(8), dimension(n_r, n_z, n_phi) :: A2can_computed, A2can_expected

        call print_test("test_compute_A2can")

        A(:,:,:,:) = 1.0d0

        call construct_zero_spline(spl_lam)
        call construct_zero_spline(spl_chi)
        call compute_A2can(A, A2can_computed)
        call destroy_splines_3d(spl_lam)
        call destroy_splines_3d(spl_chi)

        A2can_expected = A(2,:,:,:)
        if (any(abs(A2can_computed - A2can_expected) > eps)) then
            call print_fail
            print *, "A2 canonical should match A2 for identical transformation"
            error stop
        end if

        call construct_linear_spline(spl_lam)
        call construct_linear_spline(spl_chi)
        call compute_A2can(A, A2can_computed)
        call destroy_splines_3d(spl_lam)
        call destroy_splines_3d(spl_chi)

        A2can_expected = A(2,:,:,:)
        if (all(abs(A2can_computed - A2can_expected) < eps)) then
            call print_fail
            print *, "A2 canonical must not match A2 for linear transformation"
            error stop
        end if

        call print_ok
    end subroutine test_compute_A2can


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
        integer :: i_r, i_z, i_phi

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]

        do i_phi = 1, n_phi
            do i_z = 1, n_z
                do i_r = 1, n_r
                    linear(i_r, i_z, i_phi) = i_r + i_z + i_phi + 1.0d0
                end do
            end do
        end do

        call construct_splines_3d(x_min, x_max, linear, order, periodic, spl)
    end subroutine construct_linear_spline


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
