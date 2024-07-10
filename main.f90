program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: FieldType
    use magfie_tok, only: TokFieldType
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components

    implicit none

    call letter_canonical

contains

    subroutine letter_canonical
        integer, parameter :: n_r=100, n_phi=64, n_z=75
        real(dp) :: rmin, rmax, zmin, zmax
        !complex(8) :: pert

        class(FieldType), allocatable :: field_type

        ! Workaround, otherwise not initialized without perturbation field
        rmin = 75.d0
        rmax = 264.42281879194627d0
        zmin = -150.d0
        zmax = 147.38193979933115d0

        field_type = TokFieldType()

        print *, "init_canonical ..."
        call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
            [rmax, twopi, zmax], field_type)

        print *, "init_transformation ..."
        call init_transformation

        print *, "init_canonical_field_components ..."
        call init_canonical_field_components
    end subroutine letter_canonical
end program main
