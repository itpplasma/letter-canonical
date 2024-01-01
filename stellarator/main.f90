program main
    use canonical, only : init, get_transformation, spline_transformation, &
                          spl_order

    implicit none

    integer, parameter :: n_r=50, n_phi=64, n_z=75
    real(8), parameter :: rmin = 0.436d0, rmax = 2.436d0, &
                          zmin = -1.0d0, zmax = 1.0d0

    real(8), dimension(:,:,:), allocatable :: delta_phi, chi_gauge
    real(8), dimension(:,:,:,:,:,:,:), allocatable :: spl_data
    integer :: outfile_unit

    allocate(delta_phi(n_r, n_phi, n_z), chi_gauge(n_r, n_phi, n_z))
    allocate(spl_data(2,spl_order+1,spl_order+1,spl_order+1,n_r,n_z,n_phi))

    call init(n_r, n_phi, n_z)
    call get_transformation(delta_phi, chi_gauge)

    open(newunit=outfile_unit, file="delta_phi.out")
        write(outfile_unit, *) delta_phi
    close(outfile_unit)

    open(newunit=outfile_unit, file="chi_gauge.out")
        write(outfile_unit, *) chi_gauge
    close(outfile_unit)

    call spline_transformation(spl_data)

    deallocate(delta_phi, chi_gauge)
    deallocate(spl_data)


end program main
