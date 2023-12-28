module simple_cyl

    implicit none

contains

    subroutine test()
        print *, "Hello from simple_cyl"
    end subroutine test

    ! This initializes an NCSX vacuum field
    subroutine init_ncsx()
        use biotsavart, only: initialize_biotsavart

        real(8) :: extcur(10) = &
        (/ 6.52271941985300E+05, 6.51868569367400E+05, 5.37743588647300E+05, &
           2.50000000000000E-07, 2.50000000000000E-07, 2.80949750000000E+04, &
           -5.48049500000000E+04, 3.01228950000000E+04, 9.42409100000000E+04, &
           4.55138737653200E+04 /)

        call initialize_biotsavart(extcur, "c09r00")
    end subroutine init_ncsx


    subroutine eval_bfield(x, B)
        use biotsavart, only: bfield

        real(8), intent(in) :: x(3)
        real(8), intent(out) :: B(3)

        call bfield (x(1), x(2), x(3), B(1), B(2), B(3))
    end subroutine eval_bfield

    subroutine eval_afield(x, A)
        use biotsavart, only: afield

        real(8), intent(in) :: x(3)
        real(8), intent(out) :: A(3)

        call afield (x(1), x(2), x(3), A(1), A(2), A(3))
    end subroutine eval_afield


    function fieldline_direction(t, x)
        real(8), intent(in) :: t
        real(8), intent(in) :: x(3)
        real(8) :: fieldline_direction(3)

        call eval_bfield(x, fieldline_direction)

        fieldline_direction(1)=fieldline_direction(1)/(fieldline_direction(2)/x(1))
        fieldline_direction(3)=fieldline_direction(3)/(fieldline_direction(2)/x(1))
        fieldline_direction(2) = 1.0d0
    end function fieldline_direction

end module simple_cyl
