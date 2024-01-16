program main
    use canonical, only : init, get_transformation, construct_splines, spl_order

    implicit none

    integer, parameter :: n_r=50, n_z=75, n_phi=64
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    real(8), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
    real(8), dimension(:,:,:,:,:,:,:), allocatable :: spl_lam_chi
    real(8), dimension(:,:,:,:,:,:,:), allocatable :: spl_magfie  ! A2, A3, h2, h3, Bmod
    integer :: outfile_unit

    allocate(lam_phi(n_r, n_z, n_phi), chi_gauge(n_r, n_z, n_phi))
    allocate(spl_lam_chi(2,spl_order+1,spl_order+1,spl_order+1,n_r,n_z,n_phi))
    allocate(spl_magfie(5,spl_order+1,spl_order+1,spl_order+1,n_r,n_z,n_phi))

    call init(n_r, n_z, n_phi)
    call get_transformation(lam_phi, chi_gauge)

    open(newunit=outfile_unit, file="lam_phi.out")
        write(outfile_unit, *) lam_phi
    close(outfile_unit)

    open(newunit=outfile_unit, file="chi_gauge.out")
        write(outfile_unit, *) chi_gauge
    close(outfile_unit)

    spl_lam_chi(1,1,1,1,:,:,:) = lam_phi
    spl_lam_chi(2,1,1,1,:,:,:) = chi_gauge
    ! deallocate(lam_phi, chi_gauge)

    call construct_splines(spl_lam_chi)

    call test_splines

    ! TODO: Compute A, h, Bmod from splines on grid in canonical coordinates
    ! TODO: Spline A, h, Bmod on grid in canonical coordinates

    deallocate(spl_lam_chi, spl_magfie)

contains

    subroutine test_splines
        use canonical, only : eval_splines

        real(8) :: x(3)
        real(8) :: y, dy(3), d2y(6)

        x(1) = rmin
        x(2) = zmin
        x(3) = 0.d0

        call eval_splines(spl_lam_chi(1,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "y = ", y, lam_phi(1,1,1)


        call eval_splines(spl_lam_chi(2,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "y = ", y, chi_gauge(1,1,1)
    end subroutine test_splines
end program main
