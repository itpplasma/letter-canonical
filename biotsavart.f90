module biotsavart
    implicit none

    type CoilsData
        real(8), dimension(:), allocatable :: x, y, z, current
    end type CoilsData

    contains


    subroutine init_coils_data(x, y, z, current, coils)
        real(8), intent(in) :: x(:), y(:), z(:), current(:)
        type(CoilsData), intent(out) :: coils

        integer :: n_points

        n_points = size(x)
        call allocate_coils_data(coils, n_points)

        coils%x = x
        coils%y = y
        coils%z = z
        coils%current = current
    end subroutine init_coils_data


    subroutine deinit_coils_data(coils)
        type(CoilsData), intent(inout) :: coils

        call deallocate_coils_data(coils)
    end subroutine deinit_coils_data


    subroutine load_coils_file(filename, coils)
        character(len=*), intent(in) :: filename
        type(CoilsData), intent(out) :: coils

        integer :: unit
        integer :: i, n_points

        open(newunit=unit, file=filename, status="old", action="read")
        read(unit, *) n_points
        call allocate_coils_data(coils, n_points)
        do i = 1, n_points
            read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
        end do
        close(unit)

    end subroutine load_coils_file


    subroutine save_coils_file(filename, coils)
        character(len=*), intent(in) :: filename
        type(CoilsData), intent(in) :: coils

        integer :: unit
        integer :: i, n_points

        n_points = size(coils%x)
        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, *) n_points
        do i = 1, n_points
            write(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
        end do
        close(unit)
    end subroutine save_coils_file


    subroutine allocate_coils_data(coils, n_points)
        type(CoilsData), intent(out) :: coils
        integer, intent(in) :: n_points

        allocate(coils%x(n_points), coils%y(n_points), coils%z(n_points))
        allocate(coils%current(n_points))
    end subroutine allocate_coils_data


    subroutine deallocate_coils_data(coils)
        type(CoilsData), intent(inout) :: coils

        deallocate(coils%x, coils%y, coils%z, coils%current)
    end subroutine deallocate_coils_data


end module biotsavart
