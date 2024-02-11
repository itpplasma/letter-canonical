program main
    use interpolate, only : SplineData3D, construct_splines_3d, evaluate_splines_3d, destroy_splines_3d, disp
    use canonical, only : init_canonical, get_transformation, twopi, &
        compute_canonical_field_components, spl_lam, spl_chi

    implicit none
    save

    integer, parameter :: n_r=50, n_z=75, n_phi=64
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    real(8), dimension(:,:,:), allocatable :: lam_phi, chi_gauge

    integer :: outfile_unit

    allocate(lam_phi(n_r, n_z, n_phi), chi_gauge(n_r, n_z, n_phi))

    call init_canonical(n_r, n_z, n_phi)
    call get_transformation(lam_phi, chi_gauge)

    open(newunit=outfile_unit, file="lam_phi.out")
        write(outfile_unit, *) lam_phi
    close(outfile_unit)

    open(newunit=outfile_unit, file="chi_gauge.out")
        write(outfile_unit, *) chi_gauge
    close(outfile_unit)

    call construct_splines_3d( &
        x_min = [rmin, zmin, 0.d0], x_max = [rmax, zmax, twopi], y = lam_phi, &
        order = [3, 3, 3], periodic = [.False., .False., .True.], spl = spl_lam)
    deallocate(lam_phi)

    call construct_splines_3d( &
        x_min = [rmin, zmin, 0.d0], x_max = [rmax, zmax, twopi], y = chi_gauge,&
        order = [3, 3, 3], periodic = [.False., .False., .True.], spl = spl_chi)
    deallocate(chi_gauge)

    call test_splines

    ! TODO: Compute A, h, Bmod from splines on grid in canonical coordinates
    call compute_canonical_field_components
    ! TODO: Spline A, h, Bmod on grid in canonical coordinates

contains

    subroutine test_splines

        integer, parameter :: n_r_test=49, n_z_test=74, n_phi_test=63

        real(8), dimension(:,:,:), allocatable :: lam_test, chi_test
        real(8), dimension(3) :: x

        integer :: i_r, i_z, i_phi

        allocate(lam_test(n_r_test, n_z_test, n_phi_test))
        allocate(chi_test(n_r_test, n_z_test, n_phi_test))

        do i_phi = 1, n_phi_test
            do i_z = 1, n_z_test
                do i_r = 1, n_r_test
                    x = [rmin + (rmax - rmin) * (i_r - 1) / (n_r_test - 1), &
                         zmin + (zmax - zmin) * (i_z - 1) / (n_z_test - 1), &
                         twopi * (i_phi - 1) / (n_phi_test - 1)]
                    call evaluate_splines_3d( &
                        spl_lam, x, lam_test(i_r, i_z, i_phi))
                    call evaluate_splines_3d( &
                        spl_chi, x, chi_test(i_r, i_z, i_phi))
                end do
            end do
        end do

        open(newunit=outfile_unit, file="lam_spl.out")
            write(outfile_unit, *) lam_test
        close(outfile_unit)

        open(newunit=outfile_unit, file="chi_spl.out")
            write(outfile_unit, *) chi_test
        close(outfile_unit)

        deallocate(lam_test, chi_test)
    end subroutine test_splines

end program main
