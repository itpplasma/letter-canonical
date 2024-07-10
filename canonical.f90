module canonical
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: FieldType
    use interpolate, only: SplineData3D, construct_splines_3d, &
        evaluate_splines_3d, evaluate_splines_3d_der2

    implicit none
    save

    logical, parameter :: debug = .True.

    real(dp), parameter :: pi = atan(1.d0)*4.d0
    real(dp), parameter :: twopi = atan(1.d0)*8.d0

    integer, parameter :: nper = 1  ! Number of periods (TODO: make this a variable)

    real(dp) :: rmin, rmax, zmin, zmax

    real(dp) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_phi, n_z
    real(dp) :: h_r, h_phi, h_z

    ! For splines
    real(dp) :: x_min(3), x_max(3)
    integer, parameter :: order(3) = [3, 3, 3]
    logical, parameter :: periodic(3) = [.False., .True., .False.]

    ! For splining lambda (difference between canonical and cylindrical angle)
    ! and chi (gauge transformation)
    type(SplineData3D) :: spl_lam, spl_chi

    ! For splining covariant vector potential, h=B/Bmod and Bmod
    type(SplineData3D) :: spl_Bmod, spl_A1, spl_A2, spl_A3, spl_h2, spl_h3

    class(FieldType), allocatable :: magfie_type

    ! For splining psi
    real(dp) :: psi_min, psi_max
    real(dp), dimension(:,:,:), allocatable :: psi_of_x, R_of_xc
    real(dp), dimension(:), allocatable :: psi_grid
    type(SplineData3D) :: spl_R_of_xc, spl_Aphi_of_xc

contains

    subroutine init_canonical(n_r_, n_phi_, n_z_, xmin, xmax, magfie_type_)
        use magfie, only: init_magfie

        integer, intent(in) :: n_r_, n_phi_, n_z_  ! Number of grid points
        real(dp), intent(in) :: xmin(3), xmax(3)
        class(FieldType), intent(in) :: magfie_type_

        magfie_type = magfie_type_
        call magfie_type%init_magfie()

        n_r = n_r_
        n_phi = n_phi_
        n_z = n_z_
        rmin = xmin(1)
        rmax = xmax(1)
        zmin = xmin(3)
        zmax = xmax(3)

        ! Grid spacing
        h_r = (rmax-rmin)/dble(n_r-1)
        h_phi = twopi/dble(n_phi-1)
        h_z = (zmax-zmin)/dble(n_z-1)

        x_min = [rmin, 0.d0, zmin]
        x_max = [rmax, twopi, zmax]

    end subroutine init_canonical


    subroutine get_transformation(lam_phi, chi_gauge)

        real(dp), intent(inout), dimension(:,:,:) :: lam_phi, chi_gauge

        integer, parameter :: ndim=2
        real(dp), parameter :: relerr=1d-8

        real(dp), allocatable :: y(:), dy(:)
        real(dp) :: r1, r2
        integer :: i_r, i_phi, i_z, i_ctr

        i_ctr=0

        !$omp parallel private(y, dy, i_r, i_phi, i_z, r1, r2)
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

                lam_phi(1, i_phi, i_z) = 0.d0
                chi_gauge(1, i_phi, i_z) = 0.d0

                y = 0.d0

                do i_r=2,n_r
                    r1 = rmin + h_r*dble(i_r-2)
                    r2 = rmin + h_r*dble(i_r-1)

                    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_can)

                    lam_phi(i_r, i_phi, i_z) = y(1)
                    chi_gauge(i_r, i_phi, i_z) = y(2)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine get_transformation


    subroutine rh_can(r_c, y, dy)
        use magfie, only: compute_abfield

        real(dp), intent(in) :: r_c  ! plus threadprivate phi_c, z_c from module
        real(dp), dimension(2), intent(in) :: y  ! lam_phi, chi_gauge
        real(dp), dimension(2), intent(inout) :: dy
        real(dp) :: Br, Bp, Bz, Ar, Ap, Az

        call magfie_type%compute_abfield(&
            r_c, modulo(phi_c + y(1), twopi), z_c, Ar, Ap, Az, Br, Bp, Bz)

        dy(1) = -Br/Bp         ! Must still be divided by r_c for covariant Bp
        dy(2) = Ar + Ap*dy(1)  ! Here it i_r reused so that r_c cancels out
        dy(1) = dy(1) / r_c    ! Finally we divide by r_c

    end subroutine rh_can


    subroutine init_transformation
        real(dp), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
        integer :: outfile_unit

        allocate(lam_phi(n_r, n_phi, n_z), chi_gauge(n_r, n_phi, n_z))
        call get_transformation(lam_phi, chi_gauge)

        open(newunit=outfile_unit, file="lam_phi.out")
            write(outfile_unit, *) lam_phi
        close(outfile_unit)

        open(newunit=outfile_unit, file="chi_gauge.out")
            write(outfile_unit, *) chi_gauge
        close(outfile_unit)

        call construct_splines_3d(x_min, x_max, lam_phi, order, periodic, spl_lam)
        deallocate(lam_phi)

        call construct_splines_3d(x_min, x_max, chi_gauge, order, periodic, spl_chi)
        deallocate(chi_gauge)
    end subroutine init_transformation


    subroutine init_canonical_field_components
        use magfie_test, only: AMPL

        real(dp), dimension(:,:,:,:), allocatable :: xcan, xcyl, B, Acyl
        real(dp), dimension(:,:,:,:), allocatable :: hcan, Acan
        real(dp), dimension(:,:,:), allocatable :: Bmod


        allocate(xcan(3,n_r,n_phi,n_z), xcyl(3,n_r,n_phi,n_z))
        allocate(B(3,n_r,n_phi,n_z), Acyl(3,n_r,n_phi,n_z))

        call generate_regular_grid(xcan)
        call can_to_cyl(xcan, xcyl)
        call get_physical_field(xcyl, B, Acyl)

        allocate(Bmod(n_r,n_phi,n_z), hcan(2,n_r,n_phi,n_z))
        call compute_Bmod(B, Bmod)
        call compute_hcan(B, Bmod, hcan)
        call construct_splines_3d(x_min, x_max, Bmod, order, periodic, spl_Bmod)
        call construct_splines_3d(x_min, x_max, hcan(1,:,:,:), order, periodic, spl_h2)
        call construct_splines_3d(x_min, x_max, hcan(2,:,:,:), order, periodic, spl_h3)
        deallocate(hcan, Bmod)

        allocate(Acan(3,n_r,n_phi,n_z))
        call compute_Acan(Acyl, Acan)
        call construct_splines_3d(x_min, x_max, &
            Acan(1,:,:,:), [4, 5, 5], periodic, spl_A1)
        call construct_splines_3d(x_min, x_max, &
            Acan(2,:,:,:), [5, 4, 5], periodic, spl_A2)
        call construct_splines_3d(x_min, x_max, &
            Acan(3,:,:,:), [5, 5, 4], periodic, spl_A3)
        deallocate(Acan)

        deallocate(Acyl, B, xcyl, xcan)

    end subroutine init_canonical_field_components


    subroutine init_splines_with_psi
        real(dp), dimension(:,:,:), allocatable :: Aphi_of_xc
        real(dp) :: x(3)
        integer :: i_r, i_phi, i_z, debug_unit

        call init_psi_grid()
        call init_R_of_xc()

        x_min = [rmin, 0.d0, zmin]
        x_max = [rmax, twopi, zmax]

        call construct_splines_3d([psi_min, 0.0d0, zmin], [psi_max, twopi, zmax], &
            R_of_xc, order, periodic, spl_R_of_xc)

        allocate(Aphi_of_xc(n_r, n_phi, n_z))

        !$omp parallel private(i_r, i_phi, i_z, x)
        !$omp do
        do i_z = 1, n_z
            do i_phi = 1, n_phi
                do i_r = 1, n_r
                    x = get_grid_point(i_r, i_phi, i_z)
                    x(1) = R_of_xc(i_r, i_phi, i_z)
                    call evaluate_splines_3d(spl_A2, x, Aphi_of_xc(i_r, i_phi, i_z))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        call construct_splines_3d([psi_min, 0.0d0, zmin], [psi_max, twopi, zmax], &
            Aphi_of_xc, order, periodic, spl_Aphi_of_xc)

        if (debug) then
            open(newunit=debug_unit, file="Aphi_of_xc.out")
            write(debug_unit, *) Aphi_of_xc
            close(debug_unit)
        end if

    end subroutine init_splines_with_psi


    subroutine init_psi_grid
        real(dp) :: x(3)
        integer :: i_r, i_phi, i_z

        allocate(psi_of_x(n_r, n_phi, n_z), psi_grid(n_r))

        !$omp parallel private(i_r, i_phi, i_z, x)
        !$omp do
        do i_z = 1, n_z
            do i_phi = 1, n_phi
                do i_r = 1, n_r
                    x = get_grid_point(i_r, i_phi, i_z)
                    call evaluate_splines_3d(spl_A3, x, psi_of_x(i_r, i_phi, i_z))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        psi_min = minval(psi_of_x)
        psi_max = maxval(psi_of_x)

        do i_r = 1, n_r
            psi_grid(i_r) = psi_min + (psi_max - psi_min) * (i_r - 1) / (n_r - 1)
        end do
    end subroutine init_psi_grid


    subroutine init_R_of_xc
        real(dp), parameter :: tol = 1d-13
        integer, parameter :: counter_max = 32

        real(dp) :: x(3), psi, dpsi(3), d2psi(6)
        integer :: i_r, i_phi, i_z, counter

        allocate(R_of_xc(n_r, n_phi, n_z))

        !$omp parallel private(i_r, i_phi, i_z, x, psi, dpsi, d2psi, counter)
        !$omp do
        do i_z = 1, n_z
            do i_phi = 1, n_phi
                do i_r = 1, n_r
                    x = [rmin + (rmax - rmin) * (i_r - 1) / (n_r - 1), &
                         twopi * (i_phi - 1) / (n_phi - 1), &
                         zmin + (zmax - zmin) * (i_z - 1) / (n_z - 1)]
                    ! Newton iteration
                    counter = 0
                    call evaluate_splines_3d_der2(spl_A3, x, psi, dpsi, d2psi)
                    do while (abs(1.0d0 - psi_grid(i_r)/psi) > tol)
                        counter = counter + 1
                        x(1) = x(1) - (psi - psi_grid(i_r))/dpsi(1)
                        call evaluate_splines_3d_der2(spl_A3, x, psi, dpsi, d2psi)
                        if (counter >= counter_max) then
                            print *, "Newton didn't converge at", counter_max, " steps."
                            print *, counter, counter_max, x(1), psi, psi_grid(i_r)
                            error stop
                        end if
                    end do
                    R_of_xc(i_r, i_phi, i_z) = x(1)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine init_R_of_xc


    subroutine evaluate_afield_can(x, A1, dA1, d2A1, A2, dA2, d2A2, A3, dA3, d2A3)
        use magfie_test, only: AMPL

        real(dp), intent(in) :: x(3)  ! R, phi_c, Z
        real(dp), intent(out) :: A1, A2, A3, dA1(3), dA2(3), dA3(3), &
                                d2A1(6), d2A2(6), d2A3(6)

        call evaluate_splines_3d_der2(spl_A1, x, A1, dA1, d2A1)
        call evaluate_splines_3d_der2(spl_A2, x, A2, dA2, d2A2)
        call evaluate_splines_3d_der2(spl_A3, x, A3, dA3, d2A3)

    end subroutine evaluate_afield_can


    subroutine compute_Acan(Acyl, Acan)
        real(dp), dimension(:,:,:,:), intent(in) :: Acyl   ! physical components
        real(dp), dimension(:,:,:,:), intent(inout) :: Acan  ! covariant
        ! TODO: Acan only with second and third component, as the first vanishes

        integer :: i_phi, i_z, i_r
        real(dp) :: x(3)
        real(dp) :: ARcov, AZcov, Aphicov, A1can, A2can, A3can
        real(dp) :: lam, chi, dlam(3), dchi(3), dummy(6)

        do i_z=1,n_z
            do i_phi=1,n_phi
                do i_r=1,n_r
                    x = get_grid_point(i_r, i_phi, i_z)
                    call evaluate_splines_3d_der2(spl_lam, x, lam, dlam, dummy)
                    call evaluate_splines_3d_der2(spl_chi, x, chi, dchi, dummy)
                    ARcov = Acyl(1, i_r, i_phi, i_z)
                    Aphicov = Acyl(2, i_r, i_phi, i_z)*x(1)
                    AZcov = Acyl(3, i_r, i_phi, i_z)

                    A1can = ARcov + Aphicov*dlam(1) - dchi(1)
                    A2can = Aphicov*(1.0d0 + dlam(2)) - dchi(2)
                    A3can = AZcov + Aphicov*dlam(3) - dchi(3)

                    Acan(1, i_r, i_phi, i_z) = A1can
                    Acan(2, i_r, i_phi, i_z) = A2can
                    Acan(3, i_r, i_phi, i_z) = A3can
                enddo
            enddo
        enddo
    end subroutine compute_Acan


    subroutine get_physical_field(xcyl, B, A)
        use magfie_test, only: AMPL
        use magfie, only: compute_abfield

        ! Order of coordinates: R, Z, phi
        real(dp), intent(in) :: xcyl(:,:,:,:)
        real(dp), intent(inout) :: B(:,:,:,:), A(:,:,:,:)

        real(dp) :: A1, A2, A3, B1, B2, B3

        integer :: i_r, i_phi, i_z
        real(dp) :: r, z, phi

        do i_z=1,n_z
            do i_phi=1,n_phi
                do i_r=1,n_r
                    r = xcyl(1, i_r, i_phi, i_z)
                    phi = modulo(xcyl(2,i_r, i_phi, i_z), twopi)
                    z = xcyl(3, i_r, i_phi, i_z)
                    call magfie_type%compute_abfield(r, phi, z, A1, A2, A3, B1, B2, B3)
                    A(1, i_r, i_phi, i_z) = A1
                    A(2, i_r, i_phi, i_z) = A2
                    A(3 ,i_r, i_phi, i_z) = A3
                    B(1, i_r, i_phi, i_z) = B1
                    B(2, i_r, i_phi, i_z) = B2
                    B(3, i_r, i_phi, i_z) = B3
                end do
            end do
        end do

    end subroutine get_physical_field


    subroutine can_to_cyl(xcan, xcyl)
        ! TODO: write this also for vector and not just array
        real(dp), intent(in) :: xcan(:,:,:,:)
        real(dp), intent(out) :: xcyl(:,:,:,:)

        integer :: i_r, i_phi, i_z
        real(dp) :: lam

        do i_phi=1,n_phi
        do i_z=1,n_z
        do i_r=1,n_r
            call evaluate_splines_3d(spl_lam, xcan(:, i_r, i_phi, i_z), lam)
            xcyl(1, i_r, i_phi, i_z) = xcan(1, i_r, i_phi, i_z)
            xcyl(2, i_r, i_phi, i_z) = modulo(xcan(2, i_r, i_phi, i_z) + lam, twopi)
            xcyl(3, i_r, i_phi, i_z) = xcan(3, i_r, i_phi, i_z)
        enddo
        enddo
        enddo

    end subroutine can_to_cyl


    subroutine compute_Bmod(B, Bmod)
        real(dp), dimension(:,:,:,:), intent(in) :: B   ! physical components
        real(dp), dimension(:,:,:), intent(inout) :: Bmod

        integer :: i_r, i_phi, i_z

        do i_phi=1,n_phi
            do i_z=1,n_z
                do i_r=1,n_r
                    Bmod(i_r, i_phi, i_z) = sqrt( &
                        B(1, i_r, i_phi, i_z)**2 + &
                        B(2, i_r, i_phi, i_z)**2 + &
                        B(3, i_r, i_phi, i_z)**2 &
                    )
                enddo
            enddo
        enddo
    end subroutine


    subroutine compute_hcan(Bcyl, Bmod, hcan)
        real(dp), dimension(:,:,:,:), intent(in) :: Bcyl   ! physical components
        real(dp), dimension(:,:,:), intent(in) :: Bmod
        real(dp), dimension(:,:,:,:), intent(inout) :: hcan  ! covariant

        integer :: i_r, i_phi, i_z
        real(dp) :: x(3)
        real(dp) :: BRcov, BZcov, Bphicov, B1can, B2can, B3can
        real(dp) :: lam, dlam(3), dummy(6)

        do i_z=1,n_z
            do i_phi=1,n_phi
                do i_r=1,n_r
                    x = get_grid_point(i_r, i_phi, i_z)
                    call evaluate_splines_3d_der2(spl_lam, x, lam, dlam, dummy)
                    BRcov = Bcyl(1, i_r, i_phi, i_z)
                    Bphicov = Bcyl(2, i_r, i_phi, i_z)*x(1)
                    BZcov = Bcyl(3, i_r, i_phi, i_z)

                    B1can = BRcov + Bphicov*dlam(1)
                    B2can = Bphicov*(1.0d0 + dlam(2))
                    B3can = BZcov + Bphicov*dlam(3)

                    hcan(1, i_r, i_phi, i_z) = B2can/Bmod(i_r, i_phi, i_z)
                    hcan(2, i_r, i_phi, i_z) = B3can/Bmod(i_r, i_phi, i_z)
                enddo
            enddo
        enddo
    end subroutine compute_hcan


    pure subroutine cyl_to_cov(xcyl, V)
        real(dp), intent(in) :: xcyl(:,:,:,:)
        real(dp), intent(inout) :: V(:,:,:,:)

        V(2,:,:,:) = V(2,:,:,:)*xcyl(1,:,:,:)
    end subroutine cyl_to_cov


    pure subroutine generate_regular_grid(x)
        real(dp), intent(inout) :: x(:,:,:,:)

        integer :: i_r, i_phi, i_z

        do i_z=1,n_z
            do i_phi=1,n_phi
                do i_r=1,n_r
                    x(:, i_r, i_phi, i_z) = get_grid_point(i_r, i_phi, i_z)
                enddo
            enddo
        enddo

    end subroutine generate_regular_grid


    pure function get_grid_point(i_r, i_phi, i_z)
        integer, intent(in) :: i_r, i_phi, i_z
        real(dp) :: get_grid_point(3)

        get_grid_point = [ &
            rmin + h_r*dble(i_r-1), &
            h_phi*dble(i_phi-1), &
            zmin + h_z*dble(i_z-1) &
        ]
    end function get_grid_point

end module canonical
