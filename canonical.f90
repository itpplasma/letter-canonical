module canonical
    use interpolate, only : SplineData3D, evaluate_splines_3d
    use my_little_magfie, only : rmin, rmax, zmin, zmax, my_field

    implicit none

    real(8), parameter :: pi = atan(1.d0)*4.d0
    real(8), parameter :: twopi = atan(1.d0)*8.d0

    integer, parameter :: nper = 1  ! Number of periods (TODO: make thi_r a variable)

    real(8) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_z, n_phi
    real(8) :: h_r, h_z, h_phi

    ! For splining lambda (difference between canonical and cylindrical angle)
    ! and chi (gauge transformation)
    type(SplineData3D) :: spl_lam, spl_chi

    ! TODO
    ! type(SplineData3D) :: spl_A2, spl_A3, spl_h2, spl_h3, spl_Bmod

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


    subroutine compute_canonical_field_components
        real(8), dimension(:,:,:,:), allocatable :: xcan, xcyl, B, A

        allocate(xcan(3,n_r,n_z,n_phi), xcyl(3,n_r,n_z,n_phi))
        allocate(B(3,n_r,n_z,n_phi), A(3,n_r,n_z,n_phi))

        call generate_regular_grid(xcan)
        call can_to_cyl(xcan, xcyl)
        call get_physical_field(xcyl, B, A)

    end subroutine compute_canonical_field_components


    subroutine get_physical_field(xcyl, B, A)
        real(8), intent(in) :: xcyl(:,:,:,:)
        real(8), intent(out) :: B(:,:,:,:), A(:,:,:,:)

        integer :: i_r, i_z, i_phi
        real(8) :: r, z, phi

        real(8) :: dummy

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    r = xcyl(1,i_r,i_z,i_phi)
                    z = xcyl(2,i_r,i_z,i_phi)
                    phi = xcyl(3,i_r,i_z,i_phi)
                    call my_field(r, phi, z, &
                        B(1,i_r,i_z,i_phi), &
                        B(2,i_r,i_z,i_phi), &
                        B(3,i_r,i_z,i_phi), &
                        dummy, dummy, dummy, dummy, &
                        dummy, dummy, dummy, dummy, dummy, &
                        A(1,i_r,i_z,i_phi), &
                        A(2,i_r,i_z,i_phi), &
                        A(3,i_r,i_z,i_phi))
                end do
            end do
        end do

    end subroutine get_physical_field


    subroutine can_to_cyl(xcan, xcyl)
        real(8), intent(in) :: xcan(:,:,:,:)
        real(8), intent(out) :: xcyl(:,:,:,:)

        integer :: i_r, i_z, i_phi
        real(8) :: r, z, phi, lam

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    r = xcan(1,i_r,i_z,i_phi)
                    z = xcan(2,i_r,i_z,i_phi)
                    phi = xcan(3,i_r,i_z,i_phi)
                    xcyl(1,i_r,i_z,i_phi) = r
                    xcyl(2,i_r,i_z,i_phi) = z
                    call evaluate_splines_3d(spl_lam, [r, z, phi], lam)
                    xcyl(3,i_r,i_z,i_phi) = phi + lam
                enddo
            enddo
        enddo

    end subroutine can_to_cyl


    subroutine compute_A2can(A, A2can)
        real(8), dimension(:,:,:,:), intent(in) :: A
        real(8), dimension(:,:,:), intent(out) :: A2can

        integer :: i_phi, i_z, i_r

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    A2can(i_r, i_z, i_phi) = A(2, i_r, i_z, i_phi)
                enddo
            enddo
        enddo
    end subroutine compute_A2can


    subroutine generate_regular_grid(x)
        real(8), intent(out) :: x(:,:,:,:)

        integer :: i_r, i_z, i_phi

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    x(:,i_r,i_z,i_phi) = get_grid_point(i_r, i_z, i_phi)
                enddo
            enddo
        enddo

    end subroutine generate_regular_grid


    function get_grid_point(i_r, i_z, i_phi)
        integer, intent(in) :: i_r, i_z, i_phi
        real(8) :: get_grid_point(3)

        get_grid_point = [ &
            rmin + h_r*dble(i_r-1), &
            zmin + h_z*dble(i_z-1), &
            h_phi*dble(i_phi-1) &
        ]
    end function get_grid_point

end module canonical
