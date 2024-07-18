program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use magfie, only: FieldType
    use magfie_factory, only: magfie_type_from_string
    use canonical, only: twopi, init_canonical, init_transformation, &
        init_canonical_field_components
    use field_can_cyl
    use integrator

    implicit none

    integer, parameter :: n_r=100, n_phi=64, n_z=75
    integer :: nt
    character(*), parameter :: outname = "orbit.out"
    real(dp), parameter :: qe = 1d5, m = 1d0, c = 1d0, mu = 0d0 !1d-5

    class(FieldType), allocatable :: field_type
    class(field_can_cyl_t), allocatable :: field
    type(field_can_data_t) :: f
    class(symplectic_integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si

    real(dp) :: rmin, rmax, zmin, zmax
    real(dp) :: z0(4), vpar0, starttime, endtime
    real(dp), allocatable :: out(:, :)

    integer :: kt, ierr, nmax

    ! Configuration in letter_canonical.in
    character(16) :: magfie_type
    namelist /letter_canonical/ magfie_type

    nt = 8000

    ! Workaround, otherwise not initialized without perturbation field
    rmin = 75.d0
    rmax = 264.42281879194627d0
    zmin = -150.d0
    zmax = 147.38193979933115d0

    open(unit=10, file='letter_canonical.in', status='old', action='read')
    read(10, nml=letter_canonical)
    close(10)
    field_type = magfie_type_from_string(trim(magfie_type))

    print *, "init_canonical ..."
    call init_canonical(n_r, n_phi, n_z, [rmin, 0.0d0, zmin], &
        [rmax, twopi, zmax], field_type)

    print *, "init_transformation ..."
    call init_transformation

    print *, "init_canonical_field_components ..."
    call init_canonical_field_components

    ! Initial conditions
    z0(1) = 173d0   ! r
    z0(2) = 27.5d0 ! z
    z0(3) = -2.0d0*3.1415d0  ! phi
    vpar0 = 0.8d0  ! parallel velocity

    field = field_can_cyl_t()
    call field_can_init(f, mu, c*m/qe, vpar0)
    call field%evaluate(f, z0(1), z0(2), z0(3), 2)

    z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi
    call get_derivatives2(f, z0(4))
    print *, 'f%hph = ', f%hph
    print *, 'f%Aph = ', f%Aph
    print *, 'f%ro0 = ', f%ro0
    print *, 'vpar0 = ', vpar0
    print *, 'z0 = ', z0

    integ = symplectic_integrator_rk45_t(field)
    call integrator_init(si, field, f, z0, dt=1.0d-3, ntau=1, rtol=1d-8)

    allocate(out(5,nt))

    out=0d0

    out(1:4,1) = z0
    out(5,1) = f%H
    starttime = omp_get_wtime()
    nmax = nt
    do kt = 2, nt
        ierr = 0
        call integ%timestep(si, f, ierr)
        if (.not. ierr==0) then
            print *, si%z
            nmax = kt-1
            print *, 'nmax = ', nmax
            exit
        endif
        out(1:4,kt) = si%z
        out(5,kt) = f%H
    end do
    endtime = omp_get_wtime()
    print *, trim(outname), endtime-starttime

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nmax
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)

end program main
