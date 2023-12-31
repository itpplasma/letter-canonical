program main
    use canonical, only : init, get_phi_transformation

    implicit none

    integer, parameter :: n_r=50, n_phi=64, n_z=75
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    real(8), allocatable :: delta_phi(:,:,:)

    allocate(delta_phi(n_r, n_phi, n_z))

    call init
    call get_phi_transformation(rmin, rmax, zmin, zmax, delta_phi)

    write(10, *) delta_phi

end program main
