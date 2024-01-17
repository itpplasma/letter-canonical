program main
    use canonical, only : init_canonical, get_transformation, construct_splines, spl_order, twopi

    implicit none

    integer, parameter :: n_r=50, n_z=75, n_phi=64
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    real(8), dimension(:,:,:), allocatable :: lam_phi, chi_gauge, test_function
    real(8), dimension(:,:,:,:,:,:,:), allocatable :: spl_lam_chi
    real(8), dimension(:,:,:,:,:,:,:), allocatable :: spl_magfie  ! A2, A3, h2, h3, Bmod

    real(8), dimension(:), allocatable :: r, z, phi

    integer :: kr, kz, kphi, outfile_unit

    allocate(lam_phi(n_r, n_z, n_phi), chi_gauge(n_r, n_z, n_phi), test_function(n_r, n_z, n_phi))
    allocate(spl_lam_chi(3,spl_order+1,spl_order+1,spl_order+1,n_r,n_z,n_phi))
    allocate(spl_magfie(5,spl_order+1,spl_order+1,spl_order+1,n_r,n_z,n_phi))
    allocate(r(n_r), z(n_z), phi(n_phi))

    do kr = 1, n_r
        r(kr) = rmin + (rmax - rmin) * (kr - 1) / (n_r - 1)
    end do

    do kz = 1, n_z
        z(kz) = zmin + (zmax - zmin) * (kz - 1) / (n_z - 1)
    end do

    do kphi = 1, n_phi
        phi(kphi) = twopi * (kphi - 1) / n_phi
    end do

    call init_canonical(n_r, n_z, n_phi)
    call get_transformation(lam_phi, chi_gauge)

    open(newunit=outfile_unit, file="lam_phi.out")
        write(outfile_unit, *) lam_phi
    close(outfile_unit)

    open(newunit=outfile_unit, file="chi_gauge.out")
        write(outfile_unit, *) chi_gauge
    close(outfile_unit)

    do kr = 1, n_r
        do kz = 1, n_z
            do kphi = 1, n_phi
                test_function(kr, kz, kphi) = r(kr)*z(kz)*cos(phi(kphi))
            end do
        end do
    end do

    print *, "test_function(1,1,1) = ", test_function(1,1,1)

    open(newunit=outfile_unit, file="test_function.out")
        write(outfile_unit, *) test_function
    close(outfile_unit)

    spl_lam_chi(1,1,1,1,:,:,:) = lam_phi
    spl_lam_chi(2,1,1,1,:,:,:) = chi_gauge
    spl_lam_chi(3,1,1,1,:,:,:) = test_function

    ! deallocate(lam_phi, chi_gauge)

    call construct_splines(spl_lam_chi)

    call test_splines
    call test_spline_derivatives

    ! TODO: Compute A, h, Bmod from splines on grid in canonical coordinates
    ! TODO: Spline A, h, Bmod on grid in canonical coordinates

    deallocate(spl_lam_chi, spl_magfie)

contains

    subroutine test_splines
        use canonical, only : eval_splines_der2

        real(8) :: x(3)
        real(8) :: y, dy(3), d2y(6)

        integer :: kr, kz, kphi

        kr=10
        kz=20
        kphi=30

        x(1) = r(kr)
        x(2) = z(kz)
        x(3) = phi(kphi)

        call eval_splines_der2(spl_lam_chi(1,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "lam_phi = ", lam_phi(kr, kz, kphi), y


        call eval_splines_der2(spl_lam_chi(2,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "chi_gauge = ", chi_gauge(kr, kz, kphi), y

        call eval_splines_der2(spl_lam_chi(3,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "test = ", r(kr)*z(kz)*cos(phi(kphi)), y

        do kphi = 1, n_phi
            x(3) = phi(kphi)
            do kz = 1, n_z
                x(2) = z(kz)
                do kr = 1, n_r
                    x(1) = r(kr)
                    call eval_splines_der2(spl_lam_chi(1,:,:,:,:,:,:), x, y, dy, d2y)
                    lam_phi(kr, kz, kphi) = y
                    call eval_splines_der2(spl_lam_chi(2,:,:,:,:,:,:), x, y, dy, d2y)
                    chi_gauge(kr, kz, kphi) = y
                    call eval_splines_der2(spl_lam_chi(3,:,:,:,:,:,:), x, y, dy, d2y)
                    test_function(kr, kz, kphi) = y
                end do
            end do
        end do

        open(newunit=outfile_unit, file="lam_phi_splined.out")
            write(outfile_unit, *) lam_phi
        close(outfile_unit)

        open(newunit=outfile_unit, file="chi_gauge_splined.out")
            write(outfile_unit, *) chi_gauge
        close(outfile_unit)

        open(newunit=outfile_unit, file="test_function_splined.out")
            write(outfile_unit, *) test_function
        close(outfile_unit)
    end subroutine test_splines

    subroutine test_spline_derivatives
        use canonical, only : eval_splines_der2

        real(8) :: x(3)
        real(8) :: y, dy(3), d2y(6)

        integer :: kr, kz, kphi

        kr=23
        kz=32
        kphi=25

        x(1) = r(kr)
        x(2) = z(kz)
        x(3) = phi(kphi)

        call eval_splines_der2(spl_lam_chi(3,:,:,:,:,:,:), x, y, dy, d2y)

        print *, "d/dr test   = ", z(kz)*cos(phi(kphi)), dy(1)
        print *, "d/dz test   = ", r(kr)*cos(phi(kphi)), dy(2)
        print *, "d/dphi test = ", -r(kr)*z(kz)*sin(phi(kphi)), dy(3)

        print *, "d2/dr2 test   = ", 0.0, d2y(1)


    end subroutine test_spline_derivatives
end program main
