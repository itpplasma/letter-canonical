module my_little_magfie
    implicit none

contains

    ! Compute physical components of B in cylindrical coordinates x = (r,phi,z).
    subroutine eval_field_B(x, B)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3)
        real(8) :: dummy

        call field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
            dummy,dummy,dummy,dummy,dummy,dummy)
    end subroutine eval_field_B

    ! Compute physical components of B and A in cylindrical coordinates x = (r,phi,z).
    ! This routine derives from GORILLA vector_potential_rphiz (there it's covariant).
    subroutine eval_field_B_and_A(x, B, A)
        use field_eq_mod, only : rtf,btf,psif

        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3), A(3)
        real(8) :: dummy, Bmod

        call field(x(1),x(2),x(3),B(1),B(2),B(3),dummy,dummy,dummy, &
                   dummy,dummy,dummy,dummy,dummy,dummy)

        bmod=sqrt(B(1)**2+B(2)**2+B(3)**2)

        A(1)=0.d0
        A(2)=psif/x(1)
        A(3)=-rtf*btf*log(x(1))

        ! TODO: Add perturbation field
    end subroutine eval_field_B_and_A
end module my_little_magfie
