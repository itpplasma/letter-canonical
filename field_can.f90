module field_can_mod



    implicit none

    type :: FieldCan

      double precision :: A2, A3
      double precision :: h2, h3
      double precision :: Bmod

      double precision, dimension(3) :: dA2, dA3
      double precision, dimension(3) :: dh2, dh3
      double precision, dimension(3) :: dBmod

      ! second derivatives: dx1dx1, dx1dx2, dx1dx3, dx2dx2, dx2dx3, dx3dx3
      double precision, dimension(6) :: d2A2, d2A3
      double precision, dimension(6) :: d2h2, d2h3
      double precision, dimension(6) :: d2Bmod

      double precision :: H, p2, vpar
      double precision, dimension(4) :: dvpar, dH, dp2

      ! order of second derivatives:
      ! dx1dx1, dx1dx2, dx1dx3, dx2dx2, dx2dx3, dx3dx3,
      ! dpdx1, dpdx2, dpdx3, dpdp
      double precision, dimension(10) :: d2vpar, d2H, d2p2

      double precision :: mu, ro0
    end type FieldCan

    contains

    subroutine FieldCan_init(f, mu, ro0, vpar)
      type(FieldCan), intent(inout) :: f
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

    end subroutine FieldCan_init


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine get_val(f, p3)
    !
    ! computes values of H, p2 and vpar at z=(x1, x2, x3, p3)
    !
    !
      type(FieldCan), intent(inout) :: f
      double precision, intent(in) :: p3

      f%vpar = (p3 - f%A3/f%ro0)/f%h3
      f%H = f%vpar**2/2d0 + f%mu*f%Bmod  ! + qi*f%Phie
      f%p2 = f%h2*f%vpar + f%A2/f%ro0

    end subroutine get_val


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine get_derivatives(f, p3)
    !
    ! computes H, p2 and vpar at z=(x1, x2, x3, p3) and their derivatives
    !
    !
      type(FieldCan), intent(inout) :: f
      double precision, intent(in) :: p3

      call get_val(f, p3)

      f%dvpar(1:3) = -(f%dA3/f%ro0 + f%dh3*f%vpar)/f%h3
      f%dvpar(4)   = 1d0/f%h3

      f%dH(1:3) = f%vpar*f%dvpar(1:3) + f%mu*f%dBmod  ! + qi*f%dPhie
      f%dH(4)   = f%vpar/f%h3

      f%dp2(1:3) = f%dvpar(1:3)*f%h2 + f%vpar*f%dh2 + f%dA2/f%ro0
      f%dp2(4) = f%h2/f%h3

    end subroutine get_derivatives

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine get_derivatives2(f, p3)
    !
    ! computes H, p2 and vpar at z=(x1, x2, x3, p3) up to 2nd derivatives
    ! order of second derivatives:
    ! dx1dx1, dx1dx2, dx1dx3, dx2dx2, dx2dx3, dx3dx3,
    ! dpdx1, dpdx2, dpdx3, dpdp
    !
      type(FieldCan), intent(inout) :: f
      double precision, intent(in) :: p3

      call get_derivatives(f, p3)

      f%d2vpar(1:6) = -f%d2A3/f%ro0 - f%d2h3*f%vpar
      f%d2vpar(1) = f%d2vpar(1) - 2d0*f%dh3(1)*f%dvpar(1)
      f%d2vpar(2) = f%d2vpar(2) - (f%dh3(1)*f%dvpar(2) + f%dh3(2)*f%dvpar(1))
      f%d2vpar(3) = f%d2vpar(3) - (f%dh3(1)*f%dvpar(3) + f%dh3(3)*f%dvpar(1))
      f%d2vpar(4) = f%d2vpar(4) - 2d0*f%dh3(2)*f%dvpar(2)
      f%d2vpar(5) = f%d2vpar(5) - (f%dh3(2)*f%dvpar(3) + f%dh3(3)*f%dvpar(2))
      f%d2vpar(6) = f%d2vpar(6) - 2d0*f%dh3(3)*f%dvpar(3)
      f%d2vpar(1:6) = f%d2vpar(1:6)/f%h3

      f%d2H(1:6) = f%vpar*f%d2vpar(1:6) + f%mu*f%d2Bmod  ! + qi*f%d2Phie
      f%d2H(1) = f%d2H(1) + f%dvpar(1)**2
      f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
      f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)
      f%d2H(4) = f%d2H(4) + f%dvpar(2)**2
      f%d2H(5) = f%d2H(5) + f%dvpar(2)*f%dvpar(3)
      f%d2H(6) = f%d2H(6) + f%dvpar(3)**2

      f%d2p2(1:6) = f%d2vpar(1:6)*f%h2 + f%vpar*f%d2h2 + f%d2A2/f%ro0
      f%d2p2(1) = f%d2p2(1) + 2d0*f%dvpar(1)*f%dh2(1)
      f%d2p2(2) = f%d2p2(2) + f%dvpar(1)*f%dh2(2) + f%dvpar(2)*f%dh2(1)
      f%d2p2(3) = f%d2p2(3) + f%dvpar(1)*f%dh2(3) + f%dvpar(3)*f%dh2(1)
      f%d2p2(4) = f%d2p2(4) + 2d0*f%dvpar(2)*f%dh2(2)
      f%d2p2(5) = f%d2p2(5) + f%dvpar(2)*f%dh2(3) + f%dvpar(3)*f%dh2(2)
      f%d2p2(6) = f%d2p2(6) + 2d0*f%dvpar(3)*f%dh2(3)

      f%d2vpar(7:9) = -f%dh3/f%h3**2
      f%d2H(7:9) = f%dvpar(1:3)/f%h3 + f%vpar*f%d2vpar(7:9)
      f%d2p2(7:9) = f%dh2/f%h3 + f%h2*f%d2vpar(7:9)

    end subroutine get_derivatives2


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine eval_field(f, x1_c, x2_c, x3_c, mode_secders)
    !
    ! Evaluates magnetic field in canonical coordinates (x1_c, x2_c, x3_c)
    ! and stores results in variable f
    !
    ! mode_secders = 0: no second derivatives
    ! mode_secders = 1: second derivatives only in d/dx1^2
    ! mode_secders = 2: all second derivatives, including mixed
    !

      type(FieldCan), intent(inout) :: f
      double precision, intent(in) :: x1_c, x2_c, x3_c
      integer, intent(in) :: mode_secders



    end subroutine eval_field

    end module field_can_mod
