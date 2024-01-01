module my_little_magfie

    implicit none

    real(8), parameter :: rmin = 0.436d0, rmax = 2.436d0, &
                          zmin = -1.0d0, zmax = 1.0d0

contains

    subroutine test()
        print *, "Hello from Stellarator magfie"
    end subroutine test

    ! This initializes an NCSX vacuum field
    subroutine init
        use biotsavart, only: initialize_biotsavart

        real(8) :: extcur(10) = &
        (/ 6.52271941985300E+05, 6.51868569367400E+05, 5.37743588647300E+05, &
           2.50000000000000E-07, 2.50000000000000E-07, 2.80949750000000E+04, &
           -5.48049500000000E+04, 3.01228950000000E+04, 9.42409100000000E+04, &
           4.55138737653200E+04 /)

        call initialize_biotsavart(extcur, "c09r00")
    end subroutine init


    subroutine eval_field_B(x, B)
        use biotsavart, only: bfield

        real(8), intent(in) :: x(3)
        real(8), intent(out) :: B(3)

        call bfield (x(1), x(2), x(3), B(1), B(2), B(3))
    end subroutine eval_field_B

    subroutine eval_field_A(x, A)
        use biotsavart, only: afield

        real(8), intent(in) :: x(3)
        real(8), intent(out) :: A(3)

        call afield (x(1), x(2), x(3), A(1), A(2), A(3))
    end subroutine eval_field_A

    subroutine eval_field_B_and_A(x, B, A)
        real(8), intent(in) :: x(3)
        real(8), intent(inout) :: B(3), A(3)

        call eval_field_B(x, B)
        call eval_field_A(x, A)
    end subroutine eval_field_B_and_A

    subroutine my_field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az)

        real(8), intent(in) :: r, z, p
        real(8), intent(inout) :: Br,Bp,Bz,dBrdR,dBrdp,&
            dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Ar,Ap,Az

        real(8) :: x(3)
        real(8) :: B(3), A(3)

        x(1) = r
        x(2) = p
        x(3) = z

        call eval_field_B(x, B)
        call eval_field_A(x, A)

        Br = B(1)
        Bp = B(2)
        Bz = B(3)

        Ar = A(1)
        Ap = A(2)
        Az = A(3)

        dBrdR = 0.0d0
        dBrdp = 0.0d0
        dBrdZ = 0.0d0
        dBpdR = 0.0d0
        dBpdp = 0.0d0
        dBpdZ = 0.0d0
        dBzdR = 0.0d0
        dBzdp = 0.0d0
        dBzdZ = 0.0d0

    end subroutine my_field


    function fieldline_direction(t, x)
        real(8), intent(in) :: t
        real(8), intent(in) :: x(3)
        real(8) :: fieldline_direction(3)

        call eval_field_b(x, fieldline_direction)

        fieldline_direction(1) = &
            fieldline_direction(1)/(fieldline_direction(2)/x(1))
        fieldline_direction(3) = &
            fieldline_direction(3)/(fieldline_direction(2)/x(1))
        fieldline_direction(2) = 1.0d0
    end function fieldline_direction

end module my_little_magfie
