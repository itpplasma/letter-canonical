module canonical
    use my_little_magfie, only : rmin, rmax, zmin, zmax

    implicit none

    real(8), parameter :: pi = atan(1.d0)*4.d0
    integer, parameter :: spl_order = 5

    real(8) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_phi, n_z
    real(8) :: h_r, h_phi, h_z

    ! Precomputed factors for spline derivatives
    real(8), dimension(spl_order+1) :: derf1, derf2, derf3

contains

    subroutine init(n_r_, n_phi_, n_z_)
        use my_little_magfie, only : init_magfie => init

        integer, intent(in) :: n_r_, n_phi_, n_z_  ! Number of grid points
        integer :: k

        call init_magfie

        n_r = n_r_
        n_phi = n_phi_
        n_z = n_z_

        ! Grid spacing
        h_r = (rmax-rmin)/dble(n_r-1)
        h_phi = 2.d0*pi/dble(n_phi-1)
        h_z = (zmax-zmin)/dble(n_z-1)

        ! Spline derivative factors
        do k=1,spl_order+1
            derf1(k) = dble(k-1)
            derf2(k) = dble((k-1)*(k-2))
            derf3(k) = dble((k-1)*(k-2)*(k-3))
        enddo
    end subroutine init


    subroutine spline_transformation(spl_data)
        real(8), intent(inout) :: spl_data(:,:,:,:,:,:,:)

        real(8), dimension(0:spl_order,n_r)   :: splcoe_r
        real(8), dimension(0:spl_order,n_phi) :: splcoe_phi
        real(8), dimension(0:spl_order,n_z)   :: splcoe_z
        integer :: n_data
        integer :: i_r, i_p, i_z          ! Loop counters for grid
        integer :: is_z, is_phi, i_data, k  ! Loop counters for splines

        n_data = size(spl_data, 1)

        ! Spline over $\varphi$
        do i_r=1,n_r
        do i_z=1,n_z
            do i_data=1,3
                splcoe_phi(0,:)=spl_data(i_data,1,1,1,i_r,i_z,:)
                call spl_per(spl_order,n_phi,h_phi,splcoe_phi)
                do k=1,spl_order
                    spl_data(i_data,1,1,k+1,i_r,i_z,:)=splcoe_phi(k,:)
                enddo
            enddo
        enddo
        enddo

        ! Spline over $Z$
        do i_r=1,n_r
        do i_p=1,n_phi
            do is_phi=1,spl_order+1
                do i_data=1,3
                    splcoe_z(0,:)=spl_data(i_data,1,1,is_phi,i_r,:,i_p)
                    call spl_reg(spl_order,n_z,h_z,splcoe_z)
                    do k=1,spl_order
                        spl_data(i_data,1,k+1,is_phi,i_r,:,i_p)=splcoe_z(k,:)
                    enddo
                enddo
            enddo
        enddo
        enddo

        ! Spline over $R$
        do i_z=1,n_z
        do i_p=1,n_phi
            do is_z=1,spl_order+1
            do is_phi=1,spl_order+1
                do i_data=1,3
                    splcoe_r(0,:)=spl_data(i_data,1,is_z,is_phi,:,i_z,i_p)
                    call spl_reg(spl_order,n_r,h_r,splcoe_r)
                    do k=1,spl_order
                        spl_data(i_data,k+1,is_z,is_phi,:,i_z,i_p)=splcoe_r(k,:)
                    enddo
                enddo
            enddo
            enddo
        enddo
        enddo

    end subroutine spline_transformation


    subroutine get_transformation(delta_phi, chi_gauge)

        real(8), intent(inout), dimension(:,:,:) :: delta_phi, chi_gauge

        integer, parameter :: ndim=2
        real(8), parameter :: relerr=1d-10

        real(8), allocatable :: y(:), dy(:)
        real(8) :: r1, r2
        integer :: i_r, i_p, i_z, i_ctr

        i_ctr=0

        !$omp parallel private(y, dy, i_r, i_p, i_z, r1, r2)
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
            do i_p=1,n_phi
                phi_c = h_phi*dble(i_p-1)

                delta_phi(1, i_p, i_z) = 0.d0
                chi_gauge(1, i_p, i_z) = 0.d0

                y = 0.d0

                do i_r=2,n_r
                    r1 = rmin + h_r*dble(i_r-2)
                    r2 = rmin + h_r*dble(i_r-1)

                    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_ran)

                    delta_phi(i_r, i_p, i_z) = y(1)
                    chi_gauge(i_r, i_p, i_z) = y(2)
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
        dy(2) = Ar + Ap*dy(1)  ! Here it is reused so that r_c cancels out
        dy(1) = dy(1) / r_c    ! Finally we divide by r_c

    end subroutine rh_ran


end module canonical
