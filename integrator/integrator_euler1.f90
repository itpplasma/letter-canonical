module integrator_euler1
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use integrator_base, only: symplectic_integrator_t, symplectic_integrator_data_t
    use field_can, only: field_can_t, field_can_data_t, get_derivatives2

    implicit none

    type, extends(symplectic_integrator_t) :: symplectic_integrator_euler1_t
        class(field_can_t), allocatable :: field

        contains

        procedure :: timestep => orbit_timestep_sympl_euler1
    end type symplectic_integrator_euler1_t

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine orbit_timestep_sympl_euler1(self, si, f, ierr)
    !
        class(symplectic_integrator_euler1_t), intent(in) :: self
        type(symplectic_integrator_data_t), intent(inout) :: si
        type(field_can_data_t), intent(inout) :: f
        integer, intent(out) :: ierr

        integer, parameter :: n = 2
        integer, parameter :: maxit = 64

        real(dp), dimension(n) :: x, xlast
        integer :: ktau

        ! real(dp) :: x1min, x1max, x2min, x2max, fvec(2)
        ! integer :: nR, nZ, i, j

        ! nR = 100
        ! nZ = 100

        ! x1min = 163d0
        ! x1max = 164d0
        ! x2min = -28d0
        ! x2max = -27d0

        ! si%z(1) = 163.65101623987172d0
        ! si%z(2) = -27.595344619263866d0
        ! si%z(3) = 0.0d0

        ! call self%field%evaluate(f, si%z(1), si%z(2), si%z(3), 2)
        ! call get_derivatives2(f, si%z(4))
        ! x(2) = si%z(4)

        ! ! do i = 1, nR
        ! !     !do j = 1, nZ
        ! !         si%z(1) = x1min + i*(x1max-x1min)/nR
        ! !         !si%z(2) = x2min + j*(x2max-x2min)/nZ
        ! !         x(1) = si%z(1)
        ! !         x(2) = si%z(4)
        ! !         call f_sympl_euler1(si, self%field, f, n, x, fvec)
        ! !         write(6603,*) si%z, x(2), fvec
        ! !     !end do
        ! ! end do
        ! ! stop

        ierr = 0
        ktau = 0
        do while(ktau .lt. si%ntau)
            si%pthold = f%pth

            ! call self%field%evaluate(f, si%z(1), si%z(2), si%z(3), 2)
            ! call get_derivatives2(f, si%z(4))
            ! x(1) = si%z(1) - si%dt*(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
            ! x(2) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))
            x(1) = si%z(1)
            x(2) = si%z(4)

            call newton1(si, self%field, f, x, maxit, xlast, ierr)
            !call picard1(si, self%field, f, x, maxit, xlast, ierr)
            if (ierr /= 0) return

            if (x(1) < 0.0) then
                print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
                x(1) = 0.01
            end if

            si%z(1) = x(1)
            si%z(4) = x(2)

            !call extrapolate_field(f, x, xlast)
            call self%field%evaluate(f, si%z(1), si%z(2), si%z(3), 2)
            call get_derivatives2(f, si%z(4))

            !write(6602,*) si%z(1), si%z(2), si%z(3), si%z(4), si%dt*f%dH(1)/f%dpth(1)
            si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
            si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph


            ktau = ktau+1
        enddo
    end subroutine orbit_timestep_sympl_euler1


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine picard1(si, field, f, x, maxit, xlast, ierr)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer, parameter :: n = 2

        type(symplectic_integrator_data_t), intent(inout) :: si
        class(field_can_t), intent(in) :: field
        type(field_can_data_t), intent(inout) :: f
        real(dp), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(n)
        integer, intent(out) :: ierr

        integer :: kit

        ierr = 0

        do kit = 1, maxit
            call field%evaluate(f, x(1), si%z(2), si%z(3), 2)
            call get_derivatives2(f, x(2))

            x(1) = x(1) - 1d-11 *(f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2)))
            x(2) = si%z(4) - si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))/f%dpth(1)

            if(ieee_is_nan(x(1)) .or. ieee_is_nan(x(2))) then
                print *, 'picard1: nan'
                print *, si%z(1), si%z(2), si%z(3), si%z(4)
                ierr = 1
                return
            end if
            !print *, 'picard1: ', x(1), x(2)

        enddo
    end subroutine picard1

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine newton1(si, field, f, x, maxit, xlast, ierr)
    !
        integer, parameter :: n = 2

        type(symplectic_integrator_data_t), intent(inout) :: si
        class(field_can_t), intent(in) :: field
        type(field_can_data_t), intent(inout) :: f
        real(dp), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(n)
        integer, intent(out) :: ierr

        real(dp) :: fvec(n), fjac(n,n), ijac(n,n)
        real(dp) :: tolref(n)
        integer :: kit

        tolref(1) = 1d0
        tolref(2) = dabs(x(2))

        do kit = 1, maxit
            if(x(1) < 0d0) x(1) = 0.01d0

            call f_sympl_euler1(si, field, f, n, x, fvec)
            call jac_sympl_euler1(si, f, x, fjac)
            ijac(1,1) = 1d0/(fjac(1,1) - fjac(1,2)*fjac(2,1)/fjac(2,2))
            ijac(1,2) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(1,2) - fjac(2,1))
            ijac(2,1) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(2,1) - fjac(1,2))
            ijac(2,2) = 1d0/(fjac(2,2) - fjac(1,2)*fjac(2,1)/fjac(1,1))
            xlast = x
            x = x - matmul(ijac, fvec)

            ! Don't take too small values in pphi as tolerance reference
            tolref(2) = max(dabs(x(2)), tolref(2))
            tolref(2) = max(dabs(x(2)), tolref(2))

            if (all(dabs(fvec) < si%atol)) return
            if (all(dabs(x-xlast) < si%rtol*tolref)) return
        enddo
        ierr = 1
        print *, 'newton1: maximum iterations reached: ', maxit
        write(6601,*) x(1), x(2)
        write(6601,*) x-xlast
        write(6601,*) fvec
        write(6601,*) ''
        write(6601,*) fjac(1,1), fjac(1,2)
        write(6601,*) fjac(2,1), fjac(2,2)
        write(6601,*) ''
        write(6601,*) ijac(1,1), ijac(1,2)
        write(6601,*) ijac(2,1), ijac(2,2)
        write(6601,*) ''
        write(6601,*) si%z(2), si%z(3)
        write(6601,*) ''
        write(6601,*) ''
    end subroutine newton1


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine f_sympl_euler1(si, field, f, n, x, fvec)
    !
        type(symplectic_integrator_data_t), intent(inout) :: si
        class(field_can_t), intent(in) :: field
        type(field_can_data_t), intent(inout) :: f
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: fvec(n)

        call field%evaluate(f, x(1), si%z(2), si%z(3), 2)
        call get_derivatives2(f, x(2))

        fvec(1) = f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
        fvec(2) = f%dpth(1)*(x(2) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

    end subroutine f_sympl_euler1


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine jac_sympl_euler1(si, f, x, jac)
    !
        type(symplectic_integrator_data_t), intent(in) :: si
        type(field_can_data_t), intent(inout) :: f

        real(dp), intent(in)  :: x(2)
        real(dp), intent(out) :: jac(2, 2)

        jac(1,1) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)**2 &
        + si%dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
        jac(1,2) = f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
        + si%dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
        jac(2,1) = f%d2pth(1)*(x(2) - si%z(4)) &
        + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
        jac(2,2) = f%d2pth(7)*(x(2) - si%z(4)) + f%dpth(1) &
        + si%dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))

    end subroutine jac_sympl_euler1


    subroutine extrapolate_field(f, x, xlast)
        type(field_can_data_t), intent(inout) :: f
        real(dp), dimension(2), intent(in) :: x, xlast

        f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1)) + f%dpth(4)*(x(2)-xlast(2))
        f%dH(1) = f%dH(1) + f%d2H(1)*(x(1)-xlast(1)) + f%d2H(7)*(x(2)-xlast(2))
        f%dpth(1)=f%dpth(1)+f%d2pth(1)*(x(1)-xlast(1))+f%d2pth(7)*(x(2)-xlast(2))
        f%vpar = f%vpar + f%dvpar(1)*(x(1)-xlast(1)) + f%dvpar(4)*(x(2)-xlast(2))
        f%hth = f%hth + f%dhth(1)*(x(1)-xlast(1))
        f%hph = f%hph + f%dhph(1)*(x(1)-xlast(1))
    end subroutine extrapolate_field

end module integrator_euler1
