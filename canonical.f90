module canonical
    use my_little_magfie, only : rmin, rmax, zmin, zmax

    implicit none

    real(8), parameter :: pi = atan(1.d0)*4.d0
    real(8), parameter :: twopi = atan(1.d0)*8.d0

    integer, parameter :: nper = 1  ! Number of periods (TODO: make thi_r a variable)
    integer, parameter :: spl_order = 5

    real(8) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_z, n_phi
    real(8) :: h_r, h_z, h_phi

contains

    subroutine init_canonical(n_r_, n_z_, n_phi_)
        use my_little_magfie, only : init_magfie => init

        integer, intent(in) :: n_r_, n_z_, n_phi_  ! Number of grid points

        call init_magfie

        n_r = n_r_
        n_z = n_z_
        n_phi = n_phi_

        ! Grid spacing
        h_r = (rmax-rmin)/dble(n_r-1)
        h_z = (zmax-zmin)/dble(n_z-1)
        h_phi = twopi/dble(n_phi-1)

    end subroutine init_canonical


    subroutine get_transformation(delta_phi, chi_gauge)

        real(8), intent(inout), dimension(:,:,:) :: delta_phi, chi_gauge

        integer, parameter :: ndim=2
        real(8), parameter :: relerr=1d-10

        real(8), allocatable :: y(:), dy(:)
        real(8) :: r1, r2
        integer :: i_r, i_z, i_phi, i_ctr

        i_ctr=0

        !$omp parallel private(y, dy, i_r, i_z, i_phi, r1, r2)
        !$omp critical
            allocate(y(ndim),dy(ndim))
        !$omp end critical
        !$omp do
        do i_phi=1,n_phi
            !$omp critical
            i_ctr = i_ctr + 1
            print *,'integrate ODE: ',i_ctr,' of ',n_phi
            !$omp end critical
            phi_c = h_phi*dble(i_phi-1)
            do i_z=1,n_z
                z_c = zmin + h_z*dble(i_z-1)

                delta_phi(1, i_z, i_phi) = 0.d0
                chi_gauge(1, i_z, i_phi) = 0.d0

                y = 0.d0

                do i_r=2,n_r
                    r1 = rmin + h_r*dble(i_r-2)
                    r2 = rmin + h_r*dble(i_r-1)

                    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_ran)

                    delta_phi(i_r, i_z, i_phi) = y(1)
                    chi_gauge(i_r, i_z, i_phi) = y(2)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine get_transformation

    subroutine rh_ran(r_c, y, dy)
        use my_little_magfie, only : my_field

        real(8), intent(in) :: r_c  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(2), intent(in) :: y
        real(8), dimension(2), intent(inout) :: dy
        real(8) :: Br, Bp, Bz, Ar, Ap, Az
        real(8) :: dummy

        call my_field(r_c, phi_c + y(1), z_c, Br, Bp, Bz, &
            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
            Ar, Ap, Az)

        dy(1) = -Br/Bp         ! Must still be divided by r_c for covariant Bp
        dy(2) = Ar + Ap*dy(1)  ! Here it i_r reused so that r_c cancels out
        dy(1) = dy(1) / r_c    ! Finally we divide by r_c

    end subroutine rh_ran


end module canonical
