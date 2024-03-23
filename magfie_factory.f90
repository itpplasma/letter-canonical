module magfie_factory
    use magfie
    use magfie_test
    use magfie_tok

    implicit none

    contains

    function magfie_type_from_string(magfie_type) result(magfie_instance)
        character(*), intent(in) :: magfie_type
        class(FieldType), allocatable :: magfie_instance

        select case(magfie_type)
        case("test")
            magfie_instance = TestFieldType()
        case("tok")
            magfie_instance = TokFieldType()
        case default
            print *, "magfie_factory: Unknown magfie type"
            error stop
        end select
    end function magfie_type_from_string

end module magfie_factory
