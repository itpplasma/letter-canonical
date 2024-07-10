module field_can_base

  use diag_mod, only: icounter

  implicit none

  type :: field_can_data_type
    double precision :: Ath, Aph
    double precision :: hth, hph
    double precision :: Bmod

    double precision, dimension(3) :: dAth, dAph
    double precision, dimension(3) :: dhth, dhph
    double precision, dimension(3) :: dBmod

    ! second derivatives: drdr, drdth, drdph, dthdth, dthdph, dphdph
    double precision, dimension(6) :: d2Ath, d2Aph
    double precision, dimension(6) :: d2hth, d2hph
    double precision, dimension(6) :: d2Bmod

    double precision :: H, pth, vpar
    double precision, dimension(4) :: dvpar, dH, dpth

    ! order of second derivatives:
    ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
    ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
    double precision, dimension(10) :: d2vpar, d2H, d2pth

    double precision :: mu, ro0
  end type field_can_data_type

  ! abstract base class
  type, abstract :: field_can_type
    contains
    procedure(evaluate), deferred :: evaluate
  end type field_can_type

  abstract interface
    subroutine evaluate(self, f, r, th_c, ph_c, mode_secders)
      import field_can_type, field_can_data_type
      class(field_can_type), intent(in) :: self
      type(field_can_data_type), intent(inout) :: f
      double precision, intent(in) :: r, th_c, ph_c
      integer, intent(in) :: mode_secders
    end subroutine evaluate
  end interface

  contains

  subroutine FieldCan_init(f, mu, ro0, vpar, field_type)
    type(FieldCan), intent(inout) :: f
    double precision, intent(in), optional  :: mu, ro0, vpar
    integer, intent(in), optional :: field_type

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

    if (present(field_type)) then
      f%field_type = field_type
    else
      f%field_type = 0
    end if

  end subroutine FieldCan_init


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine get_val(f, pphi)
  !
  ! computes values of H, pth and vpar at z=(r, th, ph, pphi)
  !
  !
    type(FieldCan), intent(inout) :: f
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
    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: pphi

    call get_val(f, pphi)

    f%dvpar(1:3) = -(f%dAph/f%ro0 + f%dhph*f%vpar)/f%hph
    f%dvpar(4)   = 1d0/f%hph

    f%dH(1:3) = f%vpar*f%dvpar(1:3) + f%mu*f%dBmod
    f%dH(4)   = f%vpar/f%hph

    f%dpth(1:3) = f%dvpar(1:3)*f%hth + f%vpar*f%dhth + f%dAth/f%ro0

    f%dpth(4) = f%hth/f%hph

  end subroutine get_derivatives

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine get_derivatives2(f, pphi)
  !
  ! computes H, pth and vpar at z=(r, th, ph, pphi) up to 2nd derivatives
  ! order of second derivatives:
  ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
  ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
  !
    type(FieldCan), intent(inout) :: f
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


  subroutine eval_field_can_cyl(f, r, th_c, ph_c, mode_secders)
    use interpolate, only: evaluate_splines_3d_der2
    use canonical, only: spl_A2, spl_A3, spl_h2, spl_h3, spl_Bmod

    integer, parameter :: reorder(3) = [1, 3, 2]  ! dr, dph, dth -> dr, dth, dph
    integer, parameter :: reorder2(6) = [1, 3, 2, 6, 5, 4]
    ! drdr, drdph, drdth, dphdph, dphdth, dthdth ->
    ! drdr, drdth, drdph, dthdth, dthdph, dphdph

    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(8) :: x(3), a, da(3), d2a(6)

    x = [r, ph_c, th_c]

    call evaluate_splines_3d_der2(spl_A2, x, a, da, d2a)
    f%Aph = a
    f%dAph = da(reorder)
    f%d2Aph = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_A3, x, a, da, d2a)
    f%Ath = a
    f%dAth = da(reorder)
    f%d2Ath = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_h2, x, a, da, d2a)
    f%hph = a
    f%dhph = da(reorder)
    f%d2hph = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_h3, x, a, da, d2a)
    f%hth = a
    f%dhth = da(reorder)
    f%d2hth = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_Bmod, x, a, da, d2a)
    f%Bmod = a
    f%dBmod = da(reorder)
    f%d2Bmod = d2a(reorder2)

  end subroutine eval_field_can_cyl


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine eval_field(f, r, th_c, ph_c, mode_secders)
  !
  ! Evaluates magnetic field in canonical coordinates (r, th_c, ph_c)
  ! and stores results in variable f
  ! Works for A_th linear in r (toroidal flux as radial variable)
  !
  ! mode_secders = 0: no second derivatives
  ! mode_secders = 1: second derivatives only in d/dr^2
  ! mode_secders = 2: all second derivatives, including mixed
  !
  ! tested in test_magfie.f90, 2018-10-23, C. Albert <albert@alumni.tugraz.at>
  !

    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    select case (f%field_type)
      case (3)
        call eval_field_can_cyl(f, r, th_c, ph_c, mode_secders)
      case default
        print *, 'Illegal field type ', f%field_type, ' for eval_field'
        stop
    end select

  end subroutine eval_field

  end module field_can_base
