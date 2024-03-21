program main
    use interpolate, only: SplineData3D, construct_splines_3d, evaluate_splines_3d, destroy_splines_3d, disp
    use canonical, only: init_canonical, init_transformation, twopi, &
        init_canonical_field_components, spl_lam, spl_chi

    implicit none
    save

    integer, parameter :: n_r=50, n_z=75, n_phi=64
    real(8), parameter :: rmin = 75.d0, rmax = 264.42281879194627d0, &
                          zmin = -150.d0, zmax = 147.38193979933115d0

    integer :: outfile_unit

    print *, "init_canonical ..."
    call init_canonical(n_r, n_z, n_phi)

    print *, "init_transformation ..."
    call init_transformation

    print *, "init_canonical_field_components ..."
    call init_canonical_field_components

    print *, "test_splines ..."
    call test_splines

    print *, "test_integration ..."
    call test_integration


contains

    subroutine test_integration
        real(8), parameter :: tmax = 1.0d0
        integer, parameter :: nt = 10000

        real(8) :: x(3)
        integer :: i_t

        x = [0.5*(rmin+rmax), 0.5*(zmin+zmax), 0.0d0]
        do i_t = 0, nt
            call odeint_allroutines(&
                x, 3, i_t*tmax/nt, (i_t+1)*tmax/nt, 1.0d-8, Bnoncan)
            write(100, *) x
        end do

        x = [0.5*(rmin+rmax), 0.5*(zmin+zmax), 0.0d0]
        do i_t = 0, nt
            call odeint_allroutines(&
                x, 3, i_t*tmax/nt, (i_t+1)*tmax/nt, 1.0d-8, Bcan)
            write(101, *) x
        end do
    end subroutine test_integration


    subroutine Bnoncan(t, x, dx)
        use my_little_magfie, only: my_field

        real(8), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(3), intent(in) :: x
        real(8), dimension(3), intent(inout) :: dx
        real(8) :: Br, Bp, Bz, Ar, Ap, Az
        real(8) :: dummy

        call my_field(x(1), x(3), x(2), Br, Bp, Bz, &
            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
            Ar, Ap, Az)

        dx(1) = Br/Bp
        dx(2) = Bz/Bp
        dx(3) = 1.0d0
    end subroutine Bnoncan


    subroutine Bcan(t, x, dx)
        use canonical, only: spl_A2, spl_A3
        use interpolate, only: evaluate_splines_3d_der2

        real(8), intent(in) :: t  ! plus threadprivate phi_c, z_c from module
        real(8), dimension(3), intent(in) :: x
        real(8), dimension(3), intent(inout) :: dx
        real(8) :: A2, A3, dA2(3), dA3(3), d2A2(6), d2A3(6)
        real(8) :: B(3)

        call evaluate_splines_3d_der2(spl_A2, x, A2, dA2, d2A2)
        call evaluate_splines_3d_der2(spl_A3, x, A3, dA3, d2A3)

        B(1) = dA3(2) - dA2(3)
        B(2) = -dA3(1)
        B(3) = dA2(1)

        dx(1) = B(1)/B(3)
        dx(2) = B(2)/B(3)
        dx(3) = 1.0d0
    end subroutine Bcan


    subroutine test_splines

        integer, parameter :: n_r_test=49, n_z_test=74, n_phi_test=63

        real(8), dimension(:,:,:), allocatable :: lam_test, chi_test
        real(8), dimension(3) :: x

        integer :: i_r, i_z, i_phi

        allocate(lam_test(n_r_test, n_z_test, n_phi_test))
        allocate(chi_test(n_r_test, n_z_test, n_phi_test))

        do i_phi = 1, n_phi_test
            do i_z = 1, n_z_test
                do i_r = 1, n_r_test
                    x = [rmin + (rmax - rmin) * (i_r - 1) / (n_r_test - 1), &
                         zmin + (zmax - zmin) * (i_z - 1) / (n_z_test - 1), &
                         twopi * (i_phi - 1) / (n_phi_test - 1)]
                    call evaluate_splines_3d( &
                        spl_lam, x, lam_test(i_r, i_z, i_phi))
                    call evaluate_splines_3d( &
                        spl_chi, x, chi_test(i_r, i_z, i_phi))
                end do
            end do
        end do

        open(newunit=outfile_unit, file="lam_spl.out")
            write(outfile_unit, *) lam_test
        close(outfile_unit)

        open(newunit=outfile_unit, file="chi_spl.out")
            write(outfile_unit, *) chi_test
        close(outfile_unit)

        deallocate(lam_test, chi_test)
    end subroutine test_splines

end program main
