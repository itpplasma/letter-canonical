module field_can_cyl
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t, field_can_data_t

    implicit none

    type, extends(field_can_t) :: field_can_cyl_t
    contains
        procedure :: evaluate => evaluate_cyl
    end type field_can_cyl_t

    contains

    subroutine evaluate_cyl(self, f, r, th_c, ph_c, mode_secders)
        use interpolate, only: evaluate_splines_3d_der2
        use canonical, only: spl_A2, spl_A3, spl_h2, spl_h3, spl_Bmod

        class(field_can_cyl_t), intent(in) :: self
        type(field_can_data_t), intent(inout) :: f

        integer, parameter :: reorder(3) = [1, 3, 2]  ! dr, dph, dth -> dr, dth, dph
        integer, parameter :: reorder2(6) = [1, 3, 2, 6, 5, 4]
        ! drdr, drdph, drdth, dphdph, dphdth, dthdth ->
        ! drdr, drdth, drdph, dthdth, dthdph, dphdph

        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        real(dp) :: x(3), a, da(3), d2a(6)

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
    end subroutine evaluate_cyl

end module field_can_cyl
