module psi_transform
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use binsrc_sub
    use plag_coeff_sub
    implicit none

    contains

    subroutine grid_R_to_psi( n_r, n_phi, n_z, rmin, rmax, psi_inner, psi_outer, &
        psi_of_x, Aph_of_x, hph_of_x, hth_of_x, Bmod_of_x, &
        R_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc)

        integer, intent(in) :: n_r, n_phi, n_z
        real(dp), intent(in) :: rmin, rmax, psi_inner, psi_outer
        real(dp), dimension(n_r, n_phi, n_z), intent(in) :: &
          psi_of_x, Aph_of_x, hph_of_x, hth_of_x, Bmod_of_x
        real(dp), dimension(n_r, n_phi, n_z), intent(out) :: &
            R_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc

! ***** This should be an input variable:
!        integer,  intent(in) :: n_psi
! in case n_psi is different from n_r, n_r should be replaced
! with n_psi in the above declaration of output arrays
        real(dp) :: h_psi
        integer  :: n_psi
!
        integer, parameter  :: nder = 0 !(no derivatives)
! ***** This should be the input variable:
!        integer, intent(in) :: nplagr
! ***** For the moment it is a local parameter:
        integer, parameter :: nplagr = 6 !(4 points - cubic polynomial, 6: quintic)
!
        logical :: reverse
        integer :: i_R, i_psi, ibeg, iend, nshift
        integer :: i_Z, i_phi
        real(dp) :: psi_fix
        real(dp), dimension(:),   allocatable :: xp, psi_loc, R
        real(dp), dimension(:,:), allocatable :: coef
!
        allocate(xp(nplagr),coef(0:nder,nplagr),psi_loc(n_R),R(n_R))
        nshift = nplagr/2

        do i_R = 1, n_R
          R(i_R) = rmin + dble(i_R-1)*(rmax-rmin)/dble(n_R-1)
        enddo
!
! Equidistant grid in psi:
!
        n_psi = n_R   !This needs not to be the same*****
!
        i_phi = n_phi/2
        i_Z = n_Z/2
        if(psi_of_x(n_R,i_phi,i_Z).gt.psi_of_x(1,i_phi,i_Z)) then
          reverse = .false.
        else
          reverse = .true.
        endif

        h_psi=(psi_outer-psi_inner)/dble(n_psi-1)
!
        do i_phi = 1,n_phi
          do i_Z = 1,n_Z
            if(reverse) then
              psi_loc = psi_of_x(n_R:1:-1, i_phi, i_Z)
            else
              psi_loc = psi_of_x(:, i_phi, i_Z)
            endif
            do i_psi = 1,n_psi
              psi_fix = psi_inner+h_psi*dble(i_psi-1)
!
              call binsrc(psi_loc,1,n_R,psi_fix,i_R)
!
              ibeg = i_R - nshift
              iend = ibeg+nplagr-1
              if(ibeg.lt.1) then
                ibeg = 1
                iend = min(ibeg+nplagr-1,n_R)
              elseif(iend.gt.n_R) then
                iend = n_R
                ibeg = max(1, iend-nplagr+1)
              endif
!
              call plag_coeff(nplagr,nder,psi_fix,psi_loc(ibeg:iend),coef)
!
              if(reverse) then
                ibeg = n_R - ibeg +1
                iend = n_R - iend +1
                R_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*R(ibeg:iend:-1))
                Aph_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*Aph_of_x(ibeg:iend:-1, i_phi, i_Z))
                hph_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*hph_of_x(ibeg:iend:-1, i_phi, i_Z))
                hth_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*hth_of_x(ibeg:iend:-1, i_phi, i_Z))
                Bmod_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*Bmod_of_x(ibeg:iend:-1, i_phi, i_Z))
              else
                R_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*R(ibeg:iend))
                Aph_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*Aph_of_x(ibeg:iend, i_phi, i_Z))
                hph_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*hph_of_x(ibeg:iend, i_phi, i_Z))
                hth_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*hth_of_x(ibeg:iend, i_phi, i_Z))
                Bmod_of_xc(i_psi, i_phi, i_Z) = sum(coef(0,:)*Bmod_of_x(ibeg:iend, i_phi, i_Z))
              endif
            enddo
          enddo
        enddo
!
        deallocate(xp,coef,psi_loc)
    end subroutine grid_R_to_psi

end module psi_transform
