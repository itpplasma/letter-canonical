program test_field_can
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: FieldType
    use magfie_tok, only: TokFieldType
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components

    use integrator
    use field_can
    use field_can_cyl

    implicit none

    real(dp) :: z0(4), vpar0
    class(field_can_cyl_t), allocatable :: field
    type(field_can_data_t) :: f

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    real(dp) :: rmin, rmax, zmin, zmax
    !complex(8) :: pert

    class(FieldType), allocatable :: field_type

    ! Workaround, otherwise not initialized without perturbation field
    rmin = 75.d0
    rmax = 264.42281879194627d0
    zmin = -150.d0
    zmax = 147.38193979933115d0

    field_type = TokFieldType()

    print *, "init_canonical ..."
    call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
        [rmax, twopi, zmax], field_type)

    print *, "init_transformation ..."
    call init_transformation

    print *, "init_canonical_field_components ..."
    call init_canonical_field_components

    ! Initial conditions
    z0(1) = 0.1d0  ! r
    z0(2) = 0.7d0  ! theta
    z0(3) = 0.1d0  ! phi
    vpar0 = 0.8d0  ! parallel velocity

    field = field_can_cyl_t()

    call field_can_init(f, 1d-5, 1d0, vpar0)
    call field%evaluate(f, z0(1), z0(2), z0(3), 0)

    call do_test

contains

    function relerr(a, b)
        real(dp) :: relerr
        real(dp), intent(in) :: a, b
        relerr = merge(a, (a - b)/b, b == 0d0)
    end function relerr


    subroutine der2(x0, pphi, i, j)
        real(dp), intent(in) :: x0(3)
        integer, intent(in) :: i, j
        real(dp) hi, hj
        type(field_can_data_t) :: f00, f01, f10, f11
        type(field_can_data_t) :: d2fnum
        real(dp) :: pphi, x(3), dxi(3), dxj(3)
        real(dp), dimension(10) ::  d2vparnum, d2Hnum, d2pthnum
        real(dp) :: vpar00, vpar11, vpar10, vpar01, &
            H00, H11, H10, H01, pth00, pth11, pth10, pth01
        integer :: k

        hi = 1d-4
        hj = 1d-4

        dxi = 0d0
        dxj = 0d0
        dxi(i) = .5d0*hi
        dxj(j) = .5d0*hj

        x = x0 - dxi - dxj
        call field%evaluate(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f00 = f
        vpar00 = f%vpar
        H00 = f%H
        pth00 = f%pth
        x = x0 - dxi + dxj
        call field%evaluate(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f01 = f
        vpar01 = f%vpar
        H10 = f%H
        pth01 = f%pth
        x = x0 + dxi - dxj
        call field%evaluate(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f10 = f
        vpar10 = f%vpar
        H01 = f%H
        pth10 = f%pth
        x = x0 + dxi + dxj
        call field%evaluate(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f11 = f
        vpar11 = f%vpar
        H11 = f%H
        pth11 = f%pth

        call field%evaluate(f, x0(1), x0(2), x0(3), 2)
        call get_derivatives2(f, pphi)

        if (i==1 .and. j==1) k=1
        if ((i==1 .and. j==2) .or. (i==2 .and. j==1)) k=2
        if ((i==1 .and. j==3) .or. (i==3 .and. j==1)) k=3

        if (i==2 .and. j==2) k=4
        if ((i==2 .and. j==3) .or. (i==3 .and. j==2)) k=5

        if (i==3 .and. j==3) k=6

        if(i==j) then
            d2fnum%d2Ath(k) = (f11%Ath - 2d0*f%Ath + f00%Ath)/(hi*hj)
            d2fnum%d2Aph(k) = (f11%Aph - 2d0*f%Aph + f00%Aph)/(hi*hj)
            d2fnum%d2hth(k) = (f11%hth - 2d0*f%hth + f00%hth)/(hi*hj)
            d2fnum%d2hph(k) = (f11%hph - 2d0*f%hph + f00%hph)/(hi*hj)
            d2fnum%d2Bmod(k) = (f11%Bmod - 2d0*f%Bmod + f00%Bmod)/(hi*hj)
            d2vparnum(k) = (vpar11 - 2d0*f%vpar + vpar00)/(hi*hj)
            d2Hnum(k) = (H11 - 2d0*f%H + H00)/(hi*hj)
            d2pthnum(k) = (pth11 - 2d0*f%pth + pth00)/(hi*hj)
        else
            d2fnum%d2Ath(k) = (f11%Ath - f10%Ath - f01%Ath + f00%Ath)/(hi*hj)
            d2fnum%d2Aph(k) = (f11%Aph - f10%Aph - f01%Aph + f00%Aph)/(hi*hj)
            d2fnum%d2hth(k) = (f11%hth - f10%hth - f01%hth + f00%hth)/(hi*hj)
            d2fnum%d2hph(k) = (f11%hph - f10%hph - f01%hph + f00%hph)/(hi*hj)
            d2fnum%d2Bmod(k) = (f11%Bmod - f10%Bmod - f01%Bmod + f00%Bmod)/(hi*hj)
            d2vparnum(k) = (vpar11 - vpar10 - vpar01 + vpar00)/(hi*hj)
            d2Hnum(k) = (H11 - H10 - H01 + H00)/(hi*hj)
            d2pthnum(k) = (pth11 - pth10 - pth01 + pth00)/(hi*hj)
        end if

        print *, 'd2Ath (',i,j,')', f%d2Ath(k), d2fnum%d2Ath(k), relerr(d2fnum%d2Ath(k), f%d2Ath(k))
        print *, 'd2Aph (',i,j,')', f%d2Aph(k), d2fnum%d2Aph(k), relerr(d2fnum%d2Aph(k), f%d2Aph(k))
        print *, 'd2hth (',i,j,')', f%d2hth(k), d2fnum%d2hth(k), relerr(d2fnum%d2hth(k), f%d2hth(k))
        print *, 'd2hph (',i,j,')', f%d2hph(k), d2fnum%d2hph(k), relerr(d2fnum%d2hph(k), f%d2hph(k))
        print *, 'd2Bmod(',i,j,')', f%d2Bmod(k), d2fnum%d2Bmod(k), relerr(d2fnum%d2Bmod(k), f%d2Bmod(k))
        print *, 'd2vpar(',i,j,')', f%d2vpar(k), d2vparnum(k), relerr(d2vparnum(k), f%d2vpar(k))
        print *, 'd2H(',i,j,')', f%d2H(k), d2Hnum(k), relerr(d2Hnum(k), f%d2H(k))
        print *, 'd2pth(',i,j,')', f%d2pth(k), d2pthnum(k), relerr(d2pthnum(k), f%d2pth(k))
    end subroutine der2


    subroutine do_test()
        real(dp) :: dz(4)
        integer :: i, j, k
        real(dp) :: dx
        type(field_can_data_t) :: dfnum
        real(dp) :: dvparnum(4), dHnum(4), dpthnum(4)

        print *, 'f\t', 'derivative\t', 'numerical derivative\t', 'relative error'

        ! quantities to test: Ath, Aph, hth, hph, Bmod, vpar, H, pth

        do k = 1,3
            dz = 0d0
            dx = 1d-8
            dz(k) = .5d0*dx
            call field%evaluate(f, z0(1) + dz(1), z0(2) + dz(2), z0(3) + dz(3), 0)
            call get_val(f, z0(4))
            dfnum%dAth(k) = f%Ath
            dfnum%dAph(k) = f%Aph
            dfnum%dhth(k) = f%hth
            dfnum%dhph(k) = f%hph
            dfnum%dBmod(k) = f%Bmod
            dvparnum(k) = f%vpar
            dHnum(k) = f%H
            dpthnum(k) = f%pth
            call field%evaluate(f, z0(1) - dz(1), z0(2) - dz(2), z0(3) - dz(3), 0)
            call get_val(f, z0(4))
            dfnum%dAth(k) = (dfnum%dAth(k) - f%Ath)/dx
            dfnum%dAph(k) = (dfnum%dAph(k) - f%Aph)/dx
            dfnum%dhth(k) = (dfnum%dhth(k) - f%hth)/dx
            dfnum%dhph(k) = (dfnum%dhph(k) - f%hph)/dx
            dfnum%dBmod(k) = (dfnum%dBmod(k) - f%Bmod)/dx
            dvparnum(k) = (dvparnum(k) - f%vpar)/dx
            dHnum(k) = (dHnum(k) - f%H)/dx
            dpthnum(k) = (dpthnum(k) - f%pth)/dx
            call field%evaluate(f, z0(1), z0(2), z0(3), 0)
            call get_derivatives(f, z0(4))

            print *, 'dAth (',k,')', f%dAth(k), dfnum%dAth(k), relerr(dfnum%dAth(k), f%dAth(k))
            print *, 'dAph (',k,')', f%dAph(k), dfnum%dAph(k), relerr(dfnum%dAph(k), f%dAph(k))
            print *, 'dhth (',k,')', f%dhth(k), dfnum%dhth(k), relerr(dfnum%dhth(k), f%dhth(k))
            print *, 'dhph (',k,')', f%dhph(k), dfnum%dhph(k), relerr(dfnum%dhph(k), f%dhph(k))
            print *, 'dBmod(',k,')', f%dBmod(k), dfnum%dBmod(k), relerr(dfnum%dBmod(k), f%dBmod(k))
        enddo

        dx = 1d-8*z0(4)
        call get_val(f, z0(4) + .5d0*dx)
        dvparnum(4) = f%vpar
        dHnum(4) = f%H
        dpthnum(4) = f%pth
        call get_val(f, z0(4) - .5d0*dx)
        dvparnum(4) = (dvparnum(4) - f%vpar)/dx
        dHnum(4) = (dHnum(4) - f%H)/dx
        dpthnum(4) = (dpthnum(4) - f%pth)/dx
        call get_derivatives(f, z0(4))

        do k=1,3
            print *, 'dvpar(',k,')', f%dvpar(k), dvparnum(k), relerr(dvparnum(k), f%dvpar(k))
            print *, 'dH   (',k,')', f%dH(k), dHnum(k), relerr(dHnum(k), f%dH(k))
            print *, 'dpth (',k,')', f%dpth(k), dpthnum(k), relerr(dpthnum(k), f%dpth(k))
        enddo
        print *, 'dvpardpph', f%dvpar(4), dvparnum(4), relerr(dvparnum(4), f%dvpar(4))
        print *, 'dHdpph   ', f%dH(4), dHnum(4), relerr(dHnum(4), f%dH(4))
        print *, 'dpthdpph ', f%dpth(4), dpthnum(4), relerr(dpthnum(4), f%dpth(4))

        do i = 1,3
            do j = 1,3
                call der2(z0(1:3), z0(4), i, j)
            enddo
        enddo

        ! TODO: second ders in pphi and mixed
    end subroutine do_test

end program test_field_can
