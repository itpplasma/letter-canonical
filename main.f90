program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie, only: FieldType
    use magfie_tok, only: TokFieldType
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components
    use field_can_cyl
    use integrator

    implicit none

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer, parameter :: nt=100
    character(*), parameter :: outname = "euler1.out"

    class(FieldType), allocatable :: field_type
    class(field_can_cyl_t), allocatable :: field
    type(field_can_data_t) :: f
    class(symplectic_integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si

    real(dp) :: rmin, rmax, zmin, zmax
    real(dp) :: z0(4), vpar0
    real(dp), allocatable :: out(:, :)

    real(dp) :: starttime, endtime
    integer :: kt, ierr

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

    field = field_can_cyl_t()
    call field_can_init(f, 1d-5, 1d0, vpar0)

    ! Initial conditions
    z0(1) = 170d0  ! r
    z0(2) = 20d0   ! z
    z0(3) = 0.0d0  ! phi
    vpar0 = 0.8d0  ! parallel velocity

    integ = symplectic_integrator_euler1_t()
    call integrator_init(si, field, f, z0, dt=1d-3, ntau=1, rtol=1d-13)

    allocate(out(5,nt))

    out(:,1:)=0d0

    out(1:4,1) = z0
    out(5,1) = f%H
    do kt = 2, nt
        ierr = 0
        call integ%timestep(si, f, ierr)
        if (.not. ierr==0) then
            print *, si%z
            exit
        endif
        out(1:4,kt) = si%z
        out(5,kt) = f%H
    end do
    print *, outname(1:10), endtime-starttime

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)

end program main
