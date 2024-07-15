module psi_transform
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    contains

    subroutine grid_R_to_psi( n_r, n_phi, n_z, psi_of_x, &
        R_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc)

        integer, intent(in) :: n_r, n_phi, n_z
        real(dp), dimension(n_r, n_phi, n_z), intent(in) :: &
            psi_of_x(n_r, n_phi, n_z)
        real(dp), dimension(n_r, n_phi, n_z), intent(out) :: &
            R_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc

        ! TODO
    end subroutine grid_R_to_psi

end module psi_transform
