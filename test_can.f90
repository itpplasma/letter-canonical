program test_can

    ! use new_vmec_stuff_mod, only : netcdffile, multharm,ns_A,ns_s,ns_tp
    ! use parmot_mod, only : rmu,ro0
    use velo_mod,   only : isw_field_type
    use orbit_symplectic
    use field_can_mod
    use simple
    use new_vmec_stuff_mod, only: rmajor

    implicit none
    save

    real(dp) :: z0(4), vpar0
    type(Tracer) :: norb
    type(FieldCan) :: f

    integer :: npoiper2
    real(dp) :: rbig, dtau, dtaumax

    isw_field_type = 2

    ! Initial conditions
    z0(1) = 0.1d0  ! r
    z0(2) = 0.7d0  ! theta
    z0(3) = 0.1d0  ! phi
    vpar0 = 0.8d0  ! parallel velocity

    if (isw_field_type == -1) then
      call FieldCan_init(f, 1d-5, 1d0, vpar0, isw_field_type)
      call eval_field(f, z0(1), z0(2), z0(3), 0)
    else
      call init_field(norb, 'wout.nc', 5, 5, 3, 0)

      npoiper2 = 64
      rbig = rmajor*1.0d2
      dtaumax = twopi*rbig/npoiper2
      dtau = dtaumax

      call init_params(norb, 2, 4, 3.5d6, npoiper2, 1, 1d-8)  ! fusion alphas)

      ! Initial conditions
      z0(1) = 0.1d0  ! r
      z0(2) = 0.7d0  ! theta
      z0(3) = 0.1d0  ! phi
      vpar0 = 0.1d0  ! parallel velocity

      ! ro0 = mc/e*v0, different by sqrt(2) from other modules
      ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
      call FieldCan_init(f, 0d0, ro0/dsqrt(2d0), vpar0*dsqrt(2d0), isw_field_type)
      call eval_field(f, z0(1), z0(2), z0(3), 0)
      f%mu = .5d0**2*(1.d0-vpar0**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    end if

    z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi
    print *, z0(4)
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
        type(FieldCan) :: f00, f01, f10, f11
        type(FieldCan) :: d2fnum
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
        call eval_field(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f00 = f
        vpar00 = f%vpar
        H00 = f%H
        pth00 = f%pth
        x = x0 - dxi + dxj
        call eval_field(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f01 = f
        vpar01 = f%vpar
        H10 = f%H
        pth01 = f%pth
        x = x0 + dxi - dxj
        call eval_field(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f10 = f
        vpar10 = f%vpar
        H01 = f%H
        pth10 = f%pth
        x = x0 + dxi + dxj
        call eval_field(f, x(1), x(2), x(3), 0)
        call get_val(f, pphi)
        f11 = f
        vpar11 = f%vpar
        H11 = f%H
        pth11 = f%pth

        call eval_field(f, x0(1), x0(2), x0(3), 2)
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

    subroutine test_jac1(si)
      type(SymplecticIntegrator) :: si
      real(dp) :: x1(2), dx1(2), jac1(2,2), x10(2), h1(2), jac1num(2,2), fvec1(2)
      integer :: k

      h1(1) = 1d-6
      h1(2) = z0(4)*1d-6

      do k = 1,2
        dx1 = 0d0
        dx1(k) = h1(k)*0.5d0
        x10 = si%z((/1,4/)) + (/1d-4, 1d-2/)

        x1 = x10 + dx1
        call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
        jac1num(:, k) = fvec1

        x1 = x10 - dx1
        call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
        jac1num(:, k) = (jac1num(:, k) - fvec1)/h1(k)

        x1 = x10
        call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
        call jac_sympl_euler1(si, f, x1, jac1)

      end do

      print *, 'jac_sympl_euler1(1,1)', jac1(1,1), jac1num(1,1), relerr(jac1(1,1), jac1num(1,1))
      print *, 'jac_sympl_euler1(1,2)', jac1(1,2), jac1num(1,2), relerr(jac1(1,2), jac1num(1,2))
      print *, 'jac_sympl_euler1(2,1)', jac1(2,1), jac1num(2,1), relerr(jac1(2,1), jac1num(2,1))
      print *, 'jac_sympl_euler1(2,2)', jac1(2,2), jac1num(2,2), relerr(jac1(2,2), jac1num(2,2))
    end subroutine test_jac1


    subroutine test_newton(si)
      type(SymplecticIntegrator) :: si
      integer, parameter :: n = 2
      real(dp) :: x(n), fvec(n), fjac(n,n), ijac(n,n)
      integer :: k

      x = si%z((/1,4/)) + (/1d-4, 1d-2/)

      do k=1,10
        call f_sympl_euler1(si, f, n, x, fvec, 1)
        call jac_sympl_euler1(si, f, x, fjac)
        ijac(1,1) = fjac(2,2)
        ijac(1,2) = -fjac(1,2)
        ijac(2,1) = -fjac(2,1)
        ijac(2,2) = fjac(1,1)
        ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
        x = x - matmul(ijac, fvec)
      enddo

      call f_sympl_euler1(si, f, n, x, fvec, 1)

    end subroutine


    subroutine do_test()

        type(SymplecticIntegrator) :: euler1

        real(dp) :: dz(4)
        integer :: i, j, k
        real(dp) :: dx
        type(FieldCan) :: dfnum
        real(dp) :: dvparnum(4), dHnum(4), dpthnum(4)

        print *, 'f\t', 'derivative\t', 'numerical derivative\t', 'relative error'

        ! quantities to test: Ath, Aph, hth, hph, Bmod, vpar, H, pth

        do k = 1,3
            dz = 0d0
            dx = 1d-8
            dz(k) = .5d0*dx
            call eval_field(f, z0(1) + dz(1), z0(2) + dz(2), z0(3) + dz(3), 0)
            call get_val(f, z0(4))
            dfnum%dAth(k) = f%Ath
            dfnum%dAph(k) = f%Aph
            dfnum%dhth(k) = f%hth
            dfnum%dhph(k) = f%hph
            dfnum%dBmod(k) = f%Bmod
            dvparnum(k) = f%vpar
            dHnum(k) = f%H
            dpthnum(k) = f%pth
            call eval_field(f, z0(1) - dz(1), z0(2) - dz(2), z0(3) - dz(3), 0)
            call get_val(f, z0(4))
            dfnum%dAth(k) = (dfnum%dAth(k) - f%Ath)/dx
            dfnum%dAph(k) = (dfnum%dAph(k) - f%Aph)/dx
            dfnum%dhth(k) = (dfnum%dhth(k) - f%hth)/dx
            dfnum%dhph(k) = (dfnum%dhph(k) - f%hph)/dx
            dfnum%dBmod(k) = (dfnum%dBmod(k) - f%Bmod)/dx
            dvparnum(k) = (dvparnum(k) - f%vpar)/dx
            dHnum(k) = (dHnum(k) - f%H)/dx
            dpthnum(k) = (dpthnum(k) - f%pth)/dx
            call eval_field(f, z0(1), z0(2), z0(3), 0)
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

        call orbit_sympl_init(euler1, f, z0, 1.0d0, 1, 1d-12, 0, 0)
        call test_jac1(euler1)
        call test_newton(euler1)
    end subroutine do_test

    end program test_can
