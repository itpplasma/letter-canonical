program main
    use canonical, only : init, get_transformation

    implicit none

    integer, parameter :: n_r=50, n_phi=64, n_z=75
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    real(8), dimension(:,:,:), allocatable :: delta_phi, chi_gauge
    integer :: outfile_unit

    allocate(delta_phi(n_r, n_phi, n_z), chi_gauge(n_r, n_phi, n_z))

    call init(n_r, n_phi, n_z)
    call get_transformation(delta_phi, chi_gauge)

    open(newunit=outfile_unit, file="delta_phi.out")
        write(outfile_unit, *) delta_phi
    close(outfile_unit)

    open(newunit=outfile_unit, file="chi_gauge.out")
        write(outfile_unit, *) chi_gauge
    close(outfile_unit)

    deallocate(delta_phi, chi_gauge)

end program main
