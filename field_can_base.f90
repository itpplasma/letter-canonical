module field_can_base

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

end module field_can_base
