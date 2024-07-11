program test_integrator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib
    use field_can_base, only: field_can_t, field_can_data_t
    use field_can, only: field_can_init, field_can_test_t, get_derivatives2
    use integrator_base, only: symplectic_integrator_data_t
    use integrator, only: symplectic_integrator_t, symplectic_integrator_euler1_t, &
        integrator_init


    implicit none

    character(*), parameter :: outname = "euler1.out"
    real(dp), parameter :: qe = 1d0, m = 1d0, c = 1d0, mu = 1d-5
    integer :: steps_per_bounce, nbounce, nt

    integer :: ierr, kt

    real(dp) :: z0(4), vpar0, dt, taub, starttime, endtime
    real(dp), allocatable :: out(:, :)

    class(field_can_t), allocatable :: field
    type(field_can_data_t) :: f
    class(symplectic_integrator_t), allocatable :: integ
    type(symplectic_integrator_data_t) :: si


    ! Initial conditions
    z0(1) = 0.1d0  ! r
    z0(2) = 1.5d0  ! theta
    z0(3) = 0.0d0  ! phi
    vpar0 = 0.0d0  ! parallel velocity

    field = field_can_test_t()
    call field_can_init(f, mu, c*m/qe, vpar0)
    call field%evaluate(f, z0(1), z0(2), z0(3), 2)
    z0(4) = m*vpar0*f%hph + qe/c*f%Aph  ! p_phi

    call get_derivatives2(f, z0(4))
    print *, 'f%hph = ', f%hph
    print *, 'f%Aph = ', f%Aph
    print *, 'f%ro0 = ', f%ro0
    print *, 'vpar0 = ', vpar0
    print *, 'z0 = ', z0

    taub = 7800d0  ! estimated bounce time
    nbounce = 1000

    steps_per_bounce = 8
    dt = taub/steps_per_bounce
    nt = nbounce*steps_per_bounce

    print *, 'dt: ', dt
    print *, 'timesteps: ', nt

    integ = symplectic_integrator_euler1_t(field)
    call integrator_init(si, field, f, z0, dt, ntau=1, rtol=1d-12)

    allocate(out(5,nt))

    out(:,1:)=0d0

    out(1:4,1) = z0
    out(5,1) = f%H

    starttime = omp_get_wtime()
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
    endtime = omp_get_wtime()
    print *, outname(1:10), endtime-starttime

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)

end program test_integrator
