module callback
    use callback_base, only: callback_t, callback_pointer_t
    use cut_callback, only: cut_callback_t
    implicit none

    contains

    function create_callback(callback_type) result(new_callback)
        class(callback_t), allocatable :: new_callback
        character(*) :: callback_type
        select case(callback_type)
        case('cut')
            new_callback = cut_callback_t()
        case default
            error stop 'create_callback: Unknown callback type'
        end select
    end function create_callback
end module callback
