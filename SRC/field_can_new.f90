module field_can_new
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t, field_can_data_t

    implicit none

    type, extends(field_can_t) :: field_can_new_t
    contains
        procedure :: evaluate => evaluate_new
    end type field_can_new_t

    contains

    subroutine evaluate_new(self, f, r, th_c, ph_c, mode_secders)
        use interpolate, only: evaluate_splines_3d_der2
        use canonical, only: spl_Aphi_of_xc, spl_hph_of_xc, spl_hth_of_xc, spl_Bmod_of_xc

        class(field_can_new_t), intent(in) :: self
        type(field_can_data_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        integer, parameter :: reorder(3) = [1, 3, 2]  ! dr, dph, dth -> dr, dth, dph
        integer, parameter :: reorder2(6) = [1, 3, 2, 6, 5, 4]
        ! drdr, drdph, drdth, dphdph, dphdth, dthdth ->
        ! drdr, drdth, drdph, dthdth, dthdph, dphdph

        real(dp) :: x(3), a, da(3), d2a(6)

        x = [r, ph_c, th_c]

        call evaluate_splines_3d_der2(spl_Aphi_of_xc, x, a, da, d2a)
        f%Aph = a
        f%dAph = da(reorder)
        f%d2Aph = d2a(reorder2)

        f%Ath = r
        f%dAth = [1d0, 0d0, 0d0]
        f%d2Ath = 0d0

        call evaluate_splines_3d_der2(spl_hph_of_xc, x, a, da, d2a)
        f%hph = a
        f%dhph = da(reorder)
        f%d2hph = d2a(reorder2)

        call evaluate_splines_3d_der2(spl_hth_of_xc, x, a, da, d2a)
        f%hth = a
        f%dhth = da(reorder)
        f%d2hth = d2a(reorder2)

        call evaluate_splines_3d_der2(spl_Bmod_of_xc, x, a, da, d2a)
        f%Bmod = a
        f%dBmod = da(reorder)
        f%d2Bmod = d2a(reorder2)
    end subroutine evaluate_new

end module field_can_new
