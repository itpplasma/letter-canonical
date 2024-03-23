module canonical
    use magfie, only: FieldType
    use interpolate, only: SplineData3D, construct_splines_3d, &
        evaluate_splines_3d, evaluate_splines_3d_der2

    implicit none

    real(8), parameter :: pi = atan(1.d0)*4.d0
    real(8), parameter :: twopi = atan(1.d0)*8.d0

    integer, parameter :: nper = 1  ! Number of periods (TODO: make thi_r a variable)

    real(8) :: rmin, rmax, zmin, zmax

    real(8) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_z, n_phi
    real(8) :: h_r, h_z, h_phi

    ! For splines
    real(8) :: x_min(3), x_max(3)
    integer, parameter :: order(3) = [3, 3, 3]
    logical, parameter :: periodic(3) = [.False., .False., .True.]

    ! For splining lambda (difference between canonical and cylindrical angle)
    ! and chi (gauge transformation)
    type(SplineData3D) :: spl_lam, spl_chi

    ! For splining covariant vector potential, h=B/Bmod and Bmod
    type(SplineData3D) :: spl_Bmod, spl_A2, spl_A3, spl_h2, spl_h3

    class(FieldType), allocatable :: magfie_type

contains

    subroutine init_canonical(n_r_, n_z_, n_phi_, xmin, xmax, magfie_type_)
        use magfie, only: init_magfie

        integer, intent(in) :: n_r_, n_z_, n_phi_  ! Number of grid points
        real(8), intent(in) :: xmin(3), xmax(3)
        class(FieldType), intent(in) :: magfie_type_

        magfie_type = magfie_type_
        call magfie_type%init_magfie()

        n_r = n_r_
        n_z = n_z_
        n_phi = n_phi_
        rmin = xmin(1)
        rmax = xmax(1)
        zmin = xmin(2)
        zmax = xmax(2)

        ! Grid spacing
        h_r = (rmax-rmin)/dble(n_r-1)
        h_z = (zmax-zmin)/dble(n_z-1)
        h_phi = twopi/dble(n_phi-1)

        x_min = [rmin, zmin, 0.d0]
        x_max = [rmax, zmax, twopi]

    end subroutine init_canonical


    subroutine get_transformation(lam_phi, chi_gauge)

        real(8), intent(inout), dimension(:,:,:) :: lam_phi, chi_gauge

        integer, parameter :: ndim=2
        real(8), parameter :: relerr=1d-9

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

                lam_phi(1, i_z, i_phi) = 0.d0
                chi_gauge(1, i_z, i_phi) = 0.d0

                y = 0.d0

                do i_r=2,n_r
                    r1 = rmin + h_r*dble(i_r-2)
                    r2 = rmin + h_r*dble(i_r-1)

                    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_can)

                    lam_phi(i_r, i_z, i_phi) = y(1)
                    chi_gauge(i_r, i_z, i_phi) = y(2)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine get_transformation


    subroutine rh_can(r_c, y, dy)
        use magfie, only: compute_abfield

        real(8), intent(in) :: r_c  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(2), intent(in) :: y
        real(8), dimension(2), intent(inout) :: dy
        real(8) :: Br, Bp, Bz, Ar, Ap, Az

        call magfie_type%compute_abfield(&
            r_c, -phi_c + y(1), z_c, Br, Bp, Bz, Ar, Ap, Az)

        dy(1) = -Br/Bp         ! Must still be divided by r_c for covariant Bp
        dy(2) = Ar + Ap*dy(1)  ! Here it i_r reused so that r_c cancels out
        dy(1) = dy(1) / r_c    ! Finally we divide by r_c

    end subroutine rh_can


    subroutine init_transformation
        real(8), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
        integer :: outfile_unit

        allocate(lam_phi(n_r, n_z, n_phi), chi_gauge(n_r, n_z, n_phi))
        call get_transformation(lam_phi, chi_gauge)

        open(newunit=outfile_unit, file="lam_phi.out")
            write(outfile_unit, *) lam_phi
        close(outfile_unit)

        open(newunit=outfile_unit, file="chi_gauge.out")
            write(outfile_unit, *) chi_gauge
        close(outfile_unit)

        call construct_splines_3d( &
            x_min, x_max, lam_phi, order, periodic, spl_lam)
        deallocate(lam_phi)

        call construct_splines_3d( &
            x_min, x_max, chi_gauge, order, periodic, spl_chi)
        deallocate(chi_gauge)
    end subroutine init_transformation


    subroutine init_canonical_field_components
        real(8), dimension(:,:,:,:), allocatable :: xcan, xcyl, B, Acyl
        real(8), dimension(:,:,:,:), allocatable :: hcan, Acan
        real(8), dimension(:,:,:), allocatable :: Bmod


        allocate(xcan(3,n_r,n_z,n_phi), xcyl(3,n_r,n_z,n_phi))
        allocate(B(3,n_r,n_z,n_phi), Acyl(3,n_r,n_z,n_phi))

        call generate_regular_grid(xcan)
        call can_to_cyl(xcan, xcyl)
        call get_physical_field(xcyl, B, Acyl)

        allocate(Bmod(n_r,n_z,n_phi), hcan(2,n_r,n_z,n_phi))
        call compute_Bmod(B, Bmod)
        call compute_hcan(B, Bmod, hcan)
        call construct_splines_3d(x_min, x_max, Bmod, order, periodic, spl_Bmod)
        call construct_splines_3d(x_min, x_max, &
            hcan(1,:,:,:), order, periodic, spl_h2)
        call construct_splines_3d(x_min, x_max, &
            hcan(2,:,:,:), order, periodic, spl_h3)
        deallocate(hcan, Bmod)

        allocate(Acan(2,n_r,n_z,n_phi))
        call compute_Acan(Acyl, Acan)
        call construct_splines_3d(x_min, x_max, &
            Acan(1,:,:,:), order, periodic, spl_A2)
        call construct_splines_3d(x_min, x_max, &
            Acan(2,:,:,:), order, periodic, spl_A3)
        deallocate(Acan)

        deallocate(Acyl, B, xcyl, xcan)

    end subroutine init_canonical_field_components


    subroutine get_physical_field(xcyl, B, A)
        use magfie, only: compute_abfield

        ! Order of coordinates: R, Z, phi
        real(8), intent(in) :: xcyl(:,:,:,:)
        real(8), intent(inout) :: B(:,:,:,:), A(:,:,:,:)

        integer :: i_r, i_z, i_phi
        real(8) :: r, z, phi

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    r = xcyl(1,i_r,i_z,i_phi)
                    phi = xcyl(2,i_r,i_z,i_phi)
                    z = xcyl(3,i_r,i_z,i_phi)
                    call magfie_type%compute_abfield(r, phi, z, &
                        A(1,i_r,i_z,i_phi), &
                        A(2,i_r,i_z,i_phi), &
                        A(3,i_r,i_z,i_phi), &
                        B(1,i_r,i_z,i_phi), &
                        B(2,i_r,i_z,i_phi), &
                        B(3,i_r,i_z,i_phi))
                end do
            end do
        end do

    end subroutine get_physical_field


    subroutine can_to_cyl(xcan, xcyl)
        real(8), intent(in) :: xcan(:,:,:,:)
        real(8), intent(out) :: xcyl(:,:,:,:)

        integer :: i_r, i_z, i_phi
        real(8) :: r_c, z_c, phi_c, lam

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    r_c = xcan(1,i_r,i_z,i_phi)
                    z_c = xcan(2,i_r,i_z,i_phi)
                    phi_c = xcan(3,i_r,i_z,i_phi)
                    call evaluate_splines_3d(spl_lam, [r_c, z_c, phi_c], lam)
                    xcyl(1,i_r,i_z,i_phi) = r_c
                    xcyl(2,i_r,i_z,i_phi) = -phi_c + lam
                    xcyl(3,i_r,i_z,i_phi) = z_c
                enddo
            enddo
        enddo

    end subroutine can_to_cyl


    subroutine compute_Bmod(B, Bmod)
        real(8), dimension(:,:,:,:), intent(in) :: B   ! physical components
        real(8), dimension(:,:,:), intent(inout) :: Bmod

        integer :: i_r, i_z, i_phi

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    Bmod(i_r, i_z, i_phi) = sqrt( &
                        B(1, i_r, i_z, i_phi)**2 + &
                        B(2, i_r, i_z, i_phi)**2 + &
                        B(3, i_r, i_z, i_phi)**2 &
                    )
                enddo
            enddo
        enddo
    end subroutine


    subroutine compute_hcan(Bcyl, Bmod, hcan)
        real(8), dimension(:,:,:,:), intent(in) :: Bcyl   ! physical components
        real(8), dimension(:,:,:), intent(in) :: Bmod
        real(8), dimension(:,:,:,:), intent(inout) :: hcan  ! covariant

        real(8), parameter :: check_tol = 1d-2
        integer :: i_phi, i_z, i_r
        real(8) :: x(3)
        real(8) :: BRcov, BZcov, Bphicov, B1can, B2can, B3can
        real(8) :: lam, dlam(3), dummy(6)

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    x = get_grid_point(i_r, i_z, i_phi)
                    call evaluate_splines_3d_der2(spl_lam, x, lam, dlam, dummy)
                    BRcov = Bcyl(1, i_r, i_z, i_phi)
                    Bphicov = Bcyl(2, i_r, i_z, i_phi)*x(1)
                    BZcov = Bcyl(3, i_r, i_z, i_phi)

                    B1can = BRcov + Bphicov*dlam(1)
                    B2can = BZcov + Bphicov*dlam(2)
                    B3can = Bphicov*(-1.0d0 + dlam(3))

                    hcan(1, i_r, i_z, i_phi) = B2can/Bmod(i_r, i_z, i_phi)
                    hcan(2, i_r, i_z, i_phi) = B3can/Bmod(i_r, i_z, i_phi)
                enddo
            enddo
        enddo
    end subroutine compute_hcan


    subroutine compute_Acan(Acyl, Acan)
        real(8), dimension(:,:,:,:), intent(in) :: Acyl   ! physical components
        real(8), dimension(:,:,:,:), intent(inout) :: Acan  ! covariant
        ! Acan stores only second and third component, as the first vanishes

        real(8), parameter :: check_tol = 1d-2
        integer :: i_phi, i_z, i_r
        real(8) :: x(3)
        real(8) :: ARcov, AZcov, Aphicov, A1can, A2can, A3can
        real(8) :: lam, chi, dlam(3), dchi(3), dummy(6)

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    x = get_grid_point(i_r, i_z, i_phi)
                    call evaluate_splines_3d_der2(spl_lam, x, lam, dlam, dummy)
                    call evaluate_splines_3d_der2(spl_chi, x, chi, dchi, dummy)
                    ARcov = Acyl(1, i_r, i_z, i_phi)
                    Aphicov = Acyl(2, i_r, i_z, i_phi)*x(1)
                    AZcov = Acyl(3, i_r, i_z, i_phi)

                    A1can = ARcov + Aphicov*dlam(1) - dchi(1)
                    A2can = AZcov + Aphicov*dlam(2) - dchi(2)
                    A3can = Aphicov*(-1.0d0 + dlam(3)) - dchi(3)

                    Acan(1, i_r, i_z, i_phi) = A2can
                    Acan(2, i_r, i_z, i_phi) = A3can
                enddo
            enddo
        enddo
    end subroutine compute_Acan


    subroutine cyl_to_cov(xcyl, V)
        real(8), intent(in) :: xcyl(:,:,:,:)
        real(8), intent(inout) :: V(:,:,:,:)

        V(2,:,:,:) = V(2,:,:,:)*xcyl(1,:,:,:)
    end subroutine cyl_to_cov


    subroutine generate_regular_grid(x)
        real(8), intent(inout) :: x(:,:,:,:)

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
