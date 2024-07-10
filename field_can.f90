module field_can
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_can_base, only: field_can_data_type

  implicit none

contains

    subroutine field_can_init(f, mu, ro0, vpar)
        type(field_can_data_type), intent(inout) :: f
        double precision, intent(in), optional  :: mu, ro0, vpar

        if (present(mu)) then
        f%mu = mu
        else
        f%mu = 0d0
        end if

        if (present(ro0)) then
        f%ro0 = ro0
        else
        f%ro0 = 0d0
        end if

        if (present(vpar)) then
        f%vpar = vpar
        else
        f%vpar = 0d0
        end if

    end subroutine field_can_init


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine get_val(f, pphi)
    !
    ! computes values of H, pth and vpar at z=(r, th, ph, pphi)
    !
    !
        type(field_can_data_type), intent(inout) :: f
        double precision, intent(in) :: pphi

        f%vpar = (pphi - f%Aph/f%ro0)/f%hph
        f%H = f%vpar**2/2d0 + f%mu*f%Bmod
        f%pth = f%hth*f%vpar + f%Ath/f%ro0

    end subroutine get_val


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine get_derivatives(f, pphi)
    !
    ! computes H, pth and vpar at z=(r, th, ph, pphi) and their derivatives
    !
    !
        type(field_can_data_type), intent(inout) :: f
        double precision, intent(in) :: pphi

        call get_val(f, pphi)

        f%dvpar(1:3) = -(f%dAph/f%ro0 + f%dhph*f%vpar)/f%hph
        f%dvpar(4)   = 1d0/f%hph

        f%dH(1:3) = f%vpar*f%dvpar(1:3) + f%mu*f%dBmod
        f%dH(4)   = f%vpar/f%hph

        f%dpth(1:3) = f%dvpar(1:3)*f%hth + f%vpar*f%dhth + f%dAth/f%ro0

        f%dpth(4) = f%hth/f%hph

    end subroutine get_derivatives



    subroutine get_derivatives2(f, pphi)
    !
    ! computes H, pth and vpar at z=(r, th, ph, pphi) up to 2nd derivatives
    ! order of second derivatives:
    ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
    ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
    !
        type(field_can_data_type), intent(inout) :: f
        double precision, intent(in) :: pphi

        call get_derivatives(f, pphi)

        f%d2vpar(1:6) = -f%d2Aph/f%ro0 - f%d2hph*f%vpar
        f%d2vpar(1) = f%d2vpar(1) - 2d0*f%dhph(1)*f%dvpar(1)
        f%d2vpar(2) = f%d2vpar(2) - (f%dhph(1)*f%dvpar(2) + f%dhph(2)*f%dvpar(1))
        f%d2vpar(3) = f%d2vpar(3) - (f%dhph(1)*f%dvpar(3) + f%dhph(3)*f%dvpar(1))
        f%d2vpar(4) = f%d2vpar(4) - 2d0*f%dhph(2)*f%dvpar(2)
        f%d2vpar(5) = f%d2vpar(5) - (f%dhph(2)*f%dvpar(3) + f%dhph(3)*f%dvpar(2))
        f%d2vpar(6) = f%d2vpar(6) - 2d0*f%dhph(3)*f%dvpar(3)
        f%d2vpar(1:6) = f%d2vpar(1:6)/f%hph

        f%d2H(1:6) = f%vpar*f%d2vpar(1:6) + f%mu*f%d2Bmod ! + qi*d2Phie
        f%d2H(1) = f%d2H(1) + f%dvpar(1)**2
        f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
        f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)
        f%d2H(4) = f%d2H(4) + f%dvpar(2)**2
        f%d2H(5) = f%d2H(5) + f%dvpar(2)*f%dvpar(3)
        f%d2H(6) = f%d2H(6) + f%dvpar(3)**2

        f%d2pth(1:6) = f%d2vpar(1:6)*f%hth + f%vpar*f%d2hth + f%d2Ath/f%ro0
        f%d2pth(1) = f%d2pth(1) + 2d0*f%dvpar(1)*f%dhth(1)
        f%d2pth(2) = f%d2pth(2) + f%dvpar(1)*f%dhth(2) + f%dvpar(2)*f%dhth(1)
        f%d2pth(3) = f%d2pth(3) + f%dvpar(1)*f%dhth(3) + f%dvpar(3)*f%dhth(1)
        f%d2pth(4) = f%d2pth(4) + 2d0*f%dvpar(2)*f%dhth(2)
        f%d2pth(5) = f%d2pth(5) + f%dvpar(2)*f%dhth(3) + f%dvpar(3)*f%dhth(2)
        f%d2pth(6) = f%d2pth(6) + 2d0*f%dvpar(3)*f%dhth(3)

        f%d2vpar(7:9) = -f%dhph/f%hph**2
        f%d2H(7:9) = f%dvpar(1:3)/f%hph + f%vpar*f%d2vpar(7:9)
        f%d2pth(7:9) = f%dhth/f%hph + f%hth*f%d2vpar(7:9)

    end subroutine get_derivatives2

end module field_can
