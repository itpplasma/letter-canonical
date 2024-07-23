program test
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    integer, parameter :: n_r=2, n_z=3, n_phi=4
    real(dp) :: eps = 1d-10

    call setup
    call test_magfie_factory
    call test_get_grid_point
    call test_generate_regular_grid
    call test_can_to_cyl
    call test_compute_Bmod
    call test_compute_hcan

contains

    subroutine setup
        use magfie_test, only: test_field_t
        use canonical, only: init_canonical, twopi

        real(dp) :: xmin(3), xmax(3)
        xmin = [100.0d0, -150.0d0, 0.0d0]
        xmax = [200.0d0, 150.0d0, twopi]

        call print_test("Test Setup")
        call init_canonical(n_r, n_z, n_phi, xmin, xmax, test_field_t())
        call print_ok
    end subroutine setup


    subroutine test_magfie_factory
        use magfie, only: field_t
        use magfie_test, only: test_field_t
        use magfie_factory, only: magfie_type_from_string

        class(field_t), allocatable :: field_type

        call print_test("test_magfie_factory")

        field_type = magfie_type_from_string("test")

        select type (field_type)
            type is (test_field_t)
                call print_ok
            class default
                call print_fail
                print *, "unexpected field type"
                error stop
        end select
    end subroutine test_magfie_factory


    subroutine test_get_grid_point
        use canonical, only: get_grid_point, pi, rmin, rmax, zmin, zmax

        real(dp) :: x_expected(3), x_computed(3)

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
        use canonical, only: generate_regular_grid, pi, rmin, rmax, zmin, zmax

        real(dp) :: grid_computed(3, n_r, n_z, n_phi)
        real(dp) :: x_computed(3), x_expected(3)

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

        real(dp), dimension(3, n_r, n_z, n_phi) :: B
        real(dp), dimension(n_r, n_z, n_phi) :: Bmod_computed

        call print_test("test_compute_Bmod")

        B(:,:,:,:) = sqrt(1.0d0/3.0d0)

        call compute_Bmod(B, Bmod_computed)
        if (any(abs(Bmod_computed - 1.0d0) > eps)) error stop

        call print_ok
    end subroutine test_compute_Bmod


    subroutine test_compute_hcan
        use canonical, only: compute_hcan, compute_bmod, &
            generate_regular_grid, cyl_to_cov, spl_lam, Bmod, hcan
        use interpolate, only: destroy_splines_3d

        real(dp), dimension(3, n_r, n_z, n_phi) :: B, Bcov, x
        real(dp), dimension(3, n_r, n_z, n_phi) :: hcan_expected

        call print_test("test_compute_hcan")

        call generate_regular_grid(x)

        B(:,:,:,:) = sqrt(1.0d0/3.0d0)
        Bcov = B
        call cyl_to_cov(x, Bcov)

        call compute_Bmod(B, Bmod)
        hcan_expected(1,:,:,:) = Bcov(1,:,:,:)/Bmod
        hcan_expected(2,:,:,:) = Bcov(2,:,:,:)/Bmod
        hcan_expected(3,:,:,:) = Bcov(3,:,:,:)/Bmod

        call construct_zero_spline(spl_lam)
        call compute_hcan(B)
        if (any(abs(hcan - hcan_expected) > eps)) then
            call print_fail
            print *, "should match for identical transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)

        call print_ok

    end subroutine test_compute_hcan


    subroutine test_can_to_cyl
        use canonical, only: can_to_cyl, spl_lam, generate_regular_grid, twopi
        use interpolate, only: destroy_splines_3d

        real(dp), dimension(3, n_r, n_z, n_phi) :: xcan, xcyl_computed

        call generate_regular_grid(xcan)

        call print_test("test_can_to_cyl")
        call construct_zero_spline(spl_lam)
        call can_to_cyl(xcan, xcyl_computed)
        if (any(abs(modulo(xcyl_computed - xcan, twopi)) > eps)) then
            print *, xcyl_computed - xcan
            call print_fail
            print *, "should match for identical transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)

        call construct_linear_spline(spl_lam)
        call can_to_cyl(xcan, xcyl_computed)
        if (all(abs(modulo(xcyl_computed - xcan, twopi)) < eps)) then
            call print_fail
            print *, "must not match for linear transformation"
            error stop
        end if
        call destroy_splines_3d(spl_lam)

        call print_ok
    end subroutine test_can_to_cyl


    subroutine construct_zero_spline(spl)
        use canonical, only: twopi, rmin, rmax, zmin, zmax
        use interpolate, only: SplineData3D, construct_splines_3d

        type(SplineData3D), intent(out) :: spl

        real(dp) :: x_min(3), x_max(3)
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .False., .True.]

        real(dp), dimension(n_r, n_z, n_phi) :: zeros

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]
        zeros(:, :, :) = 0.0d0

        call construct_splines_3d(x_min, x_max, zeros, order, periodic, spl)
    end subroutine construct_zero_spline


    subroutine construct_linear_spline(spl)
        use canonical, only: twopi, rmin, rmax, zmin, zmax
        use interpolate, only: SplineData3D, construct_splines_3d

        type(SplineData3D), intent(out) :: spl

        real(dp) :: x_min(3), x_max(3)
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .False., .True.]

        real(dp), dimension(n_r, n_z, n_phi) :: linear

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]

        call fill_linear(linear)

        call construct_splines_3d(x_min, x_max, linear, order, periodic, spl)
    end subroutine construct_linear_spline


    subroutine fill_linear(x)
        real(dp), intent(inout) :: x(n_r, n_z, n_phi)

        integer :: i_r, i_z, i_phi

        do i_phi = 1, n_phi
            do i_z = 1, n_z
                do i_r = 1, n_r
                    x(i_r, i_z, i_phi) = &
                        0.23d0*i_r + 0.12d0*i_z + 0.123d0*i_phi + 1.1234d0
                end do
            end do
        end do
    end subroutine fill_linear

end program test
