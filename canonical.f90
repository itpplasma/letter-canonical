module canonical
    implicit none

    real(8), parameter :: pi=atan(1.d0)*4.d0

contains

    subroutine init
        use my_little_magfie, only : init_magfie => init

        call init_magfie
    end subroutine init

    subroutine get_phi_transformation(rmin, rmax, zmin, zmax, delta_phi)
        real(8), intent(in) :: rmin, rmax, zmin, zmax
        real(8), intent(inout) :: delta_phi(:, :, :)

        integer, parameter :: ndim=1
        real(8), parameter :: relerr=1d-10

        integer :: n_r, n_phi, n_z
        integer :: i_r, i_phi, i_z, i_ctr
        real(8) :: h_r, h_phi, h_z
        real(8) :: phi_c, z_c

        real(8), allocatable :: y(:), dy(:)

        real(8) :: r1, r2

        i_ctr=0

        n_r = size(delta_phi, 1)
        n_phi = size(delta_phi, 2)
        n_z = size(delta_phi, 3)

        h_r = (rmax-rmin)/dble(n_r-1)
        h_phi = 2.d0*pi/dble(n_phi-1)
        h_z = (zmax-zmin)/dble(n_z-1)

        !$omp parallel private(y, dy, i_r, i_phi, i_z, r1, r2, r, dG_dt, dG_dp)
        !$omp critical
            allocate(y(ndim),dy(ndim))
        !$omp end critical
        !$omp do
        do i_z=1,n_z
            !$omp critical
            i_ctr = i_ctr + 1
            print *,'integrate ODE: ',i_ctr,' of ',n_z
            !$omp end critical
            z_c = zmin + h_z*dble(i_z-1)
            do i_phi=1,n_phi
                phi_c = h_phi*dble(i_phi-1)

                delta_phi(1, i_z, i_phi)=0.d0

                y(1) = 0.d0

                do i_r=2,n_r
                    r1 = rmin + h_r*dble(i_r-2)
                    r2 = rmin + h_r*dble(i_r-1)

                    call odeint_allroutines(y, ndim, r1, r2, relerr, derivs)

                    delta_phi(i_r,i_z,i_phi)=y(1)
                enddo
            enddo
        enddo
        !$omp end do

        contains

            subroutine derivs(r_c, f, df)
                real(8) :: r_c
                real(8), dimension(2) :: f, df

                call rhs_canonical(r_c, phi_c, z_c, f, df)

            end subroutine derivs

    end subroutine get_phi_transformation

    subroutine rhs_canonical(r_c, phi_c, z_c, y, dy)
        use my_little_magfie, only : my_field

        real(8) :: r_c, phi_c, z_c, dummy
        real(8) :: Br, Bp, Bz, Ar, Ap, Az
        real(8), dimension(2) :: y, dy

        call my_field(r_c, phi_c + y(1), z_c, Br, Bp, Bz, &
            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
            Ar, Ap, Az)

        dy(1) = -Br/Bp         ! Must still be divided by r_c for covariant Bp
        dy(2) = Ar + Ap*dy(1)  ! Here it is reused so that r_c cancels out
        dy(1) = dy(1) / r_c    ! Finally we divide by r_c

    end subroutine rhs_canonical


end module canonical
