program test_biotsavart
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    character(len=*), parameter :: test_coils_file = "coils.test"
    character(len=*), parameter :: real_coils_file = &
        "/proj/plasma/data/W7X/COILS/coils.w7x"

    call test_load_coils_file

    contains

    subroutine test_load_coils_file
        use biotsavart, only: CoilsData, load_coils_file, deinit_coils_data

        type(CoilsData) :: coils

        call print_test("load_coils_file")

        call create_test_coils_file
        call load_coils_file(test_coils_file, coils)
        call remove_test_coils_file

        if (size(coils%x) /= 2) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if

        call deinit_coils_data(coils)

        call print_ok
    end subroutine test_load_coils_file


    subroutine init_test_coils_data(coils)
        use biotsavart, only: CoilsData, init_coils_data

        type(CoilsData), intent(out) :: coils

        real(8), parameter :: large_distance = 1.0e6
        real(8), dimension(2) :: x, y, z, current

        x = [0.0d0, 0.0d0]
        y = [0.0d0, 0.0d0]
        z = [-large_distance, large_distance]
        current = [1.0d0, 0.0d0]

        call init_coils_data(x, y, z, current, coils)
    end subroutine init_test_coils_data


    subroutine create_test_coils_file
        use biotsavart, only: CoilsData, init_coils_data, save_coils_file

        type(CoilsData) :: coils

        call init_test_coils_data(coils)
        call save_coils_file(test_coils_file, coils)
    end subroutine create_test_coils_file


    subroutine remove_test_coils_file
        call system("rm -f " // test_coils_file)
    end subroutine remove_test_coils_file

end program test_biotsavart
