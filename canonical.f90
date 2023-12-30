module canonical
    implicit none

contains

    subroutine init
        use my_little_magfie, only : init_magfie => init

        call init_magfie
    end subroutine init

end module canonical
