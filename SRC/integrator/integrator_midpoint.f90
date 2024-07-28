module integrator_midpoint
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use integrator_base, only: symplectic_integrator_t, symplectic_integrator_data_t
    use field_can, only: field_can_t, field_can_data_t, get_derivatives2

    implicit none



    type, extends(symplectic_integrator_t) :: symplectic_integrator_midpoint_t
        class(field_can_t), allocatable :: field

        contains

        procedure :: timestep => orbit_timestep_midpoint
    end type symplectic_integrator_midpoint_t

    contains


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_midpoint(self, si, f, ierr)
!
  class(symplectic_integrator_midpoint_t), intent(in) :: self
  type(symplectic_integrator_data_t), intent(inout) :: si
  type(field_can_data_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 5
  integer, parameter :: maxit = 32

  double precision, dimension(n) :: x, xlast
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1:4) = si%z
    x(5) = si%z(1)

    call newton_midpoint(si, self%field, f, x, si%atol, si%rtol, maxit, xlast)

    si%z = x(1:4)

    f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1) + x(5) - xlast(5)) &  ! d/dr
                  + f%dpth(2)*(x(2)-xlast(2)) &  ! d/dth
                  + f%dpth(3)*(x(3)-xlast(3)) &  ! d/dph
                  + f%dpth(4)*(x(4)-xlast(4))    ! d/dpph

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_midpoint


subroutine newton_midpoint(si, field, f, x, atol, rtol, maxit, xlast)
  type(symplectic_integrator_data_t), intent(inout) :: si
  class(field_can_t), intent(in) :: field
  type(field_can_data_t), intent(inout) :: f
  type(field_can_data_t) :: fmid

  integer, parameter :: n = 5
  integer :: kit

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info

  double precision :: xabs(n), tolref(n), fabs(n)


  real(dp), parameter :: twopi = atan(1.0d0)*8.0d0

  tolref(1) = 1d7
  tolref(2) = 1d2
  tolref(3) = 1d3
  tolref(4) = 1d3
  tolref(5) = 1d7

  do kit = 1, maxit
    call f_midpoint_part1(si, field, f, n, x, fvec)
    call jac_midpoint_part1(si, f, x, fjac)
    fmid = f
    call f_midpoint_part2(si, field, f, n, x, fvec)
    call jac_midpoint_part2(si, f, fmid, x, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_midpoint: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x(1), x(2), x(3), x(4), x(5), xabs, fvec
end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_midpoint_part1(si, field, f, n, x, fvec)
  !

  type(symplectic_integrator_data_t), intent(inout) :: si
  class(field_can_t), intent(in) :: field
  type(field_can_data_t), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)

    ! evaluate at midpoint
    call field%evaluate(f, x(5), 0.5d0*(x(2) + si%z(2)), 0.5d0*(x(3) + si%z(3)), 2)
    call get_derivatives2(f, 0.5d0*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) + 0.5d0*si%dt*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))

  end subroutine f_midpoint_part1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_midpoint_part2(si, field, f, n, x, fvec)
  !
  type(symplectic_integrator_data_t), intent(inout) :: si
  class(field_can_t), intent(in) :: field
  type(field_can_data_t), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)

    double precision :: dpthmid, pthdotbar

    ! save evaluation from midpoint
    dpthmid = f%dpth(1)
    pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)

    ! evaluate at endpoint
    call field%evaluate(f, x(1), x(2), x(3), 2)
    call get_derivatives2(f, x(4))
    fvec(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

  end subroutine f_midpoint_part2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_midpoint_part1(si, f, x, jac)
  !
  type(symplectic_integrator_data_t), intent(inout) :: si
  type(field_can_data_t), intent(inout) :: f
    double precision, intent(in)  :: x(5)
    double precision, intent(out) :: jac(5, 5)

    jac(2,1) = 0d0
    jac(2,5) = f%d2pth(1)*(x(2) - si%z(2)) - si%dt*f%d2H(1)
    jac(2,2:3) = 0.5d0*(f%d2pth(2:3)*(x(2) - si%z(2)) - si%dt*f%d2H(2:3))
    jac(2,2) = jac(2,2) + f%dpth(1)
    jac(2,4) = 0.5d0*(f%d2pth(7)*(x(2) - si%z(2)) - si%dt*f%d2H(7))

    jac(3,1) = 0d0
    jac(3,5) = (f%d2pth(1)*f%hph + f%dpth(1)*f%dhph(1))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(1)*f%vpar + f%dpth(1)*f%dvpar(1) &
      - f%d2H(1)*f%hth - f%dH(1)*f%dhth(1))
    jac(3,2:3) = 0.5d0*((f%d2pth(2:3)*f%hph + f%dpth(1)*f%dhph(2:3))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(2:3)*f%vpar + f%dpth(1)*f%dvpar(2:3) &
      - f%d2H(2:3)*f%hth - f%dH(1)*f%dhth(2:3)))
    jac(3,3) = jac(3,3) + f%dpth(1)*f%hph
    jac(3,4) = 0.5d0*(f%d2pth(7)*f%hph*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(7)*f%vpar + f%dpth(1)*f%dvpar(4) - f%d2H(7)*f%hth))

    jac(4,1) = 0d0
    jac(4,5) = f%d2pth(1)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
      - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
    jac(4,2) = 0.5d0*(f%d2pth(2)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(5)*f%dpth(1) + f%dH(3)*f%d2pth(2) &
      - f%d2H(2)*f%dpth(3) - f%dH(1)*f%d2pth(5)))
    jac(4,3) = 0.5d0*(f%d2pth(3)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(6)*f%dpth(1) + f%dH(3)*f%d2pth(3) &
      - f%d2H(3)*f%dpth(3) - f%dH(1)*f%d2pth(6)))
    jac(4,4) = 0.5d0*(f%d2pth(7)*(x(4) - si%z(4)) &
      + si%dt*(f%dH(3)*f%d2pth(7) + f%d2H(9)*f%dpth(1)&
      - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9)))
    jac(4,4) = jac(4,4) + f%dpth(1)

    jac(5,1) = 0d0
    jac(5,5) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(1)&
      + si%dt/2.0d0*(f%d2pth(1)*f%dH(2) + f%dpth(1)*f%d2H(2) &
      - f%d2pth(2)*f%dH(1) - f%dpth(2)*f%d2H(1))
    jac(5,2) = 0.5d0*(f%d2pth(2)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(2) &
      + si%dt/2.0d0*(f%d2pth(2)*f%dH(2) + f%dpth(1)*f%d2H(4) &
      - f%d2pth(4)*f%dH(1) - f%dpth(2)*f%d2H(2)))
    jac(5,3) = 0.5d0*(f%d2pth(3)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(3) &
      + si%dt/2.0d0*(f%d2pth(3)*f%dH(2) + f%dpth(1)*f%d2H(5) &
      - f%d2pth(5)*f%dH(1) - f%dpth(2)*f%d2H(3)))
    jac(5,4) = 0.5d0*(f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
      + si%dt/2.0d0*(f%d2pth(7)*f%dH(2) + f%dpth(1)*f%d2H(8) &
      - f%d2pth(8)*f%dH(1) - f%dpth(2)*f%d2H(7)))

end subroutine jac_midpoint_part1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_midpoint_part2(si, f, fmid, x, jac)
  !
  type(symplectic_integrator_data_t), intent(inout) :: si
  type(field_can_data_t), intent(inout) :: f
  type(field_can_data_t), intent(inout) :: fmid
    double precision, intent(in)  :: x(5)
    double precision, intent(out) :: jac(5, 5)

    ! fmid%dpth(1)*(f%pth - si%pthold) + si%dt*(fmid%dpth(1)*fmid%dH(2)-fmid%dpth(2)*fmid%dH(1))

    jac(1,1) = fmid%dpth(1)*f%dpth(1)
    jac(1,2) = 0.5d0*(fmid%d2pth(2)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(2)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(4) &
      - fmid%dpth(2)*fmid%d2H(2) - fmid%d2pth(4)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(2)
    jac(1,3) = 0.5d0*(fmid%d2pth(3)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(3)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(5) &
      - fmid%dpth(2)*fmid%d2H(3) - fmid%d2pth(5)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(3)
    jac(1,4) = 0.5d0*(fmid%d2pth(7)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(7)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(8) &
        - fmid%dpth(2)*fmid%d2H(7) - fmid%d2pth(8)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(4)
    jac(1,5) = fmid%d2pth(1)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(1)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(2) &
        - fmid%dpth(2)*fmid%d2H(1) - fmid%d2pth(2)*fmid%dH(1))

end subroutine jac_midpoint_part2



end module integrator_midpoint
