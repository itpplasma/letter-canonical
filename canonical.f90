module canonical
    use my_little_magfie, only : rmin, rmax, zmin, zmax

    implicit none

    real(8), parameter :: pi = atan(1.d0)*4.d0
    real(8), parameter :: twopi = atan(1.d0)*2.d0

    integer, parameter :: nper = 1  ! Number of periods (TODO: make thi_r a variable)
    integer, parameter :: spl_order = 5

    real(8) :: phi_c, z_c  ! Temporary variables for the ODE integration
    !$omp threadprivate(phi_c, z_c)

    ! Number and spacing of grid points
    integer :: n_r, n_z, n_phi
    real(8) :: h_r, h_z, h_phi

    ! Precomputed factors for spline derivatives
    real(8), dimension(spl_order+1) :: derf1, derf2, derf3

contains

    subroutine init(n_r_, n_z_, n_phi_)
        use my_little_magfie, only : init_magfie => init

        integer, intent(in) :: n_r_, n_z_, n_phi_  ! Number of grid points
        integer :: k

        call init_magfie

        n_r = n_r_
        n_z = n_z_
        n_phi = n_phi_

        ! Grid spacing
        h_r = (rmax-rmin)/dble(n_r-1)
        h_z = (zmax-zmin)/dble(n_z-1)
        h_phi = twopi/dble(n_phi-1)

        ! Spline derivative factors
        do k=1,spl_order+1
            derf1(k) = dble(k-1)
            derf2(k) = dble((k-1)*(k-2))
            derf3(k) = dble((k-1)*(k-2)*(k-3))
        enddo
    end subroutine init


    subroutine construct_splines(spl_data)
        real(8), intent(inout) :: spl_data(:,:,:,:,:,:,:)

        real(8), dimension(0:spl_order,n_r)   :: splcoe_r
        real(8), dimension(0:spl_order,n_z)   :: splcoe_z
        real(8), dimension(0:spl_order,n_phi) :: splcoe_phi

        integer :: n_data
        integer :: i_r, i_z, i_phi            ! Loop counters for grid
        integer :: i_r_z, i_r_phi, i_data, k  ! Loop counters for splines

        n_data = size(spl_data, 1)

        ! Spline over $\varphi$
        do i_r=1,n_r
        do i_z=1,n_z
            do i_data=1,n_data
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
        do i_phi=1,n_phi
            do i_r_phi=1,spl_order+1
                do i_data=1,n_data
                    splcoe_z(0,:)=spl_data(i_data,1,1,i_r_phi,i_r,:,i_phi)
                    call spl_reg(spl_order,n_z,h_z,splcoe_z)
                    do k=1,spl_order
                        spl_data(i_data,1,k+1,i_r_phi,i_r,:,i_phi)=splcoe_z(k,:)
                    enddo
                enddo
            enddo
        enddo
        enddo

        ! Spline over $R$
        do i_z=1,n_z
        do i_phi=1,n_phi
            do i_r_z=1,spl_order+1
            do i_r_phi=1,spl_order+1
                do i_data=1,n_data
                    splcoe_r(0,:)=spl_data(i_data,1,i_r_z,i_r_phi,:,i_z,i_phi)
                    call spl_reg(spl_order,n_r,h_r,splcoe_r)
                    do k=1,spl_order
                        spl_data(i_data,k+1,i_r_z,i_r_phi,:,i_z,i_phi)=splcoe_r(k,:)
                    enddo
                enddo
            enddo
            enddo
        enddo
        enddo

    end subroutine construct_splines


    subroutine eval_splines(spl_data, x, y, dy, d2y)
        real(8), intent(in) :: spl_data(:,:,:,:,:,:)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: y, dy(3), d2y(6)

        real(8) :: dr, dz, dphi
        integer :: nsp1                ! Spline order plus one
        integer :: i_r, i_z, i_phi, k  ! Counters

        real(8) :: qua,dqua_dr,dqua_dt,dqua_dp
        real(8) :: d2qua_dr2,d2qua_drdt,d2qua_drdp,d2qua_dt2,d2qua_dtdp,d2qua_dp2

        real(8), dimension(spl_order+1)        :: sp_all,dsp_all_ds,dsp_all_dt
        real(8), dimension(spl_order+1)        :: d2sp_all_ds2,d2sp_all_dsdt,d2sp_all_dt2
        real(8), dimension(spl_order+1, spl_order+1) :: stp_all, dstp_all_ds,d2stp_all_ds2


        dr=x(1)/h_r
        i_r=max(0,min(n_r-1,int(dr)))
        dr=(dr-dble(i_r))*h_r
        i_r=i_r+1

        dz=x(2)/h_z
        i_z=max(0,min(n_z-1,int(dz)))
        dz=(dz-dble(i_z))*h_z
        i_z=i_z+1

        dphi=modulo(x(3),twopi/dble(nper))/h_phi
        i_phi=max(0,min(n_phi-1,int(h_phi)))
        dphi=(h_phi-dble(i_phi))*h_phi
        i_phi=i_phi+1

        nsp1 = spl_order+1

        !
        ! Begin interpolation over $x1$
        !
        stp_all(1:nsp1,1:nsp1)=spl_data(nsp1,:,:,i_r,i_z,i_phi)
        dstp_all_ds(1:nsp1,1:nsp1)=stp_all(1:nsp1,1:nsp1)*derf1(nsp1)
        d2stp_all_ds2(1:nsp1,1:nsp1)=stp_all(1:nsp1,1:nsp1)*derf2(nsp1)
        !
        do k=spl_order,3,-1
            stp_all(1:nsp1,1:nsp1)=spl_data(k,:,:,i_r,i_z,i_phi) &
                + dr*stp_all(1:nsp1,1:nsp1)
            dstp_all_ds(1:nsp1,1:nsp1)=spl_data(k,:,:,i_r,i_z,i_phi)*derf1(k) &
                + dr*dstp_all_ds(1:nsp1,1:nsp1)
            d2stp_all_ds2(1:nsp1,1:nsp1)=spl_data(k,:,:,i_r,i_z,i_phi)*derf2(k) &
                + dr*d2stp_all_ds2(1:nsp1,1:nsp1)
        enddo
        !
        stp_all(1:nsp1,1:nsp1)=spl_data(1,:,:,i_r,i_z,i_phi) &
             + dr*(spl_data(2,:,:,i_r,i_z,i_phi)+dr*stp_all(1:nsp1,1:nsp1))
        dstp_all_ds(1:nsp1,1:nsp1)=spl_data(2,:,:,i_r,i_z,i_phi) &
            + dr*dstp_all_ds(1:nsp1,1:nsp1)
        !
        ! End interpolation over $x1$
        !-------------------------------
        ! Begin interpolation over $x2$
        !
        sp_all(1:nsp1)=stp_all(nsp1,1:nsp1)
        dsp_all_ds(1:nsp1)=dstp_all_ds(nsp1,1:nsp1)
        d2sp_all_ds2(1:nsp1)=d2stp_all_ds2(nsp1,1:nsp1)
        dsp_all_dt(1:nsp1)=sp_all(1:nsp1)*derf1(nsp1)
        d2sp_all_dsdt(1:nsp1)=dsp_all_ds(1:nsp1)*derf1(nsp1)
        d2sp_all_dt2(1:nsp1)=sp_all(1:nsp1)*derf2(nsp1)

        do k=spl_order,3,-1
            sp_all(1:nsp1)=stp_all(k,1:nsp1)+dz*sp_all(1:nsp1)
            dsp_all_ds(1:nsp1)=dstp_all_ds(k,1:nsp1)+dz*dsp_all_ds(1:nsp1)
            d2sp_all_ds2(1:nsp1)=d2stp_all_ds2(k,1:nsp1)+dz*d2sp_all_ds2(1:nsp1)
            dsp_all_dt(1:nsp1)=stp_all(k,1:nsp1)*derf1(k)+dz*dsp_all_dt(1:nsp1)
            d2sp_all_dsdt(1:nsp1)=dstp_all_ds(k,1:nsp1)*derf1(k) &
                + dz*d2sp_all_dsdt(1:nsp1)
            d2sp_all_dt2(1:nsp1)=stp_all(k,1:nsp1)*derf2(k)+dz*d2sp_all_dt2(1:nsp1)
        enddo

        sp_all(1:nsp1)=stp_all(1,1:nsp1) + dz*(stp_all(2,1:nsp1)+dz*sp_all(1:nsp1))
        dsp_all_ds(1:nsp1)=dstp_all_ds(1,1:nsp1) &
                            + dz*(dstp_all_ds(2,1:nsp1)+dz*dsp_all_ds(1:nsp1))
        d2sp_all_ds2(1:nsp1)=d2stp_all_ds2(1,1:nsp1) &
                            + dz*(d2stp_all_ds2(2,1:nsp1)+dz*d2sp_all_ds2(1:nsp1))
        dsp_all_dt(1:nsp1)=stp_all(2,1:nsp1)+dz*dsp_all_dt(1:nsp1)
        d2sp_all_dsdt(1:nsp1)=dstp_all_ds(2,1:nsp1)+dz*d2sp_all_dsdt(1:nsp1)
        !
        ! End interpolation over $x2$
        !--------------------------------
        ! Begin interpolation over $x3$
        !
        qua=sp_all(nsp1)
        dqua_dr=dsp_all_ds(nsp1)
        dqua_dt=dsp_all_dt(nsp1)
        dqua_dp=qua*derf1(nsp1)

        d2qua_dr2=d2sp_all_ds2(nsp1)
        d2qua_drdt=d2sp_all_dsdt(nsp1)
        d2qua_drdp=dqua_dr*derf1(nsp1)
        d2qua_dt2=d2sp_all_dt2(nsp1)
        d2qua_dtdp=dqua_dt*derf1(nsp1)
        d2qua_dp2=qua*derf2(nsp1)

        do k=spl_order,3,-1
            qua=sp_all(k)+dphi*qua
            dqua_dr=dsp_all_ds(k)+dphi*dqua_dr
            dqua_dt=dsp_all_dt(k)+dphi*dqua_dt
            dqua_dp=sp_all(k)*derf1(k)+dphi*dqua_dp

            d2qua_dr2=d2sp_all_ds2(k)+dphi*d2qua_dr2
            d2qua_drdt=d2sp_all_dsdt(k)+dphi*d2qua_drdt
            d2qua_drdp=dsp_all_ds(k)*derf1(k)+dphi*d2qua_drdp
            d2qua_dt2=d2sp_all_dt2(k)+dphi*d2qua_dt2
            d2qua_dtdp=dsp_all_dt(k)*derf1(k)+dphi*d2qua_dtdp
            d2qua_dp2=sp_all(k)*derf2(k)+dphi*d2qua_dp2
        enddo

        qua=sp_all(1)+dphi*(sp_all(2)+dphi*qua)
        dqua_dr=dsp_all_ds(1)+dphi*(dsp_all_ds(2)+dphi*dqua_dr)
        dqua_dt=dsp_all_dt(1)+dphi*(dsp_all_dt(2)+dphi*dqua_dt)

        d2qua_dr2=d2sp_all_ds2(1)+dphi*(d2sp_all_ds2(2)+dphi*d2qua_dr2)
        d2qua_drdt=d2sp_all_dsdt(1)+dphi*(d2sp_all_dsdt(2)+dphi*d2qua_drdt)
        d2qua_dt2=d2sp_all_dt2(1)+dphi*(d2sp_all_dt2(2)+dphi*d2qua_dt2)

        dqua_dp=sp_all(2)+dphi*dqua_dp
        d2qua_drdp=dsp_all_ds(2)+dphi*d2qua_drdp
        d2qua_dtdp=dsp_all_dt(2)+dphi*d2qua_dtdp
        !
        ! End interpolation over $x3$
        !

        y=qua

        dy(1)=dqua_dr
        dy(2)=dqua_dt
        dy(3)=dqua_dp

        d2y(1)=d2qua_dr2
        d2y(2)=d2qua_drdt
        d2y(3)=d2qua_drdp
        d2y(4)=d2qua_dt2
        d2y(5)=d2qua_dtdp
        d2y(6)=d2qua_dp2
    end subroutine eval_splines


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
