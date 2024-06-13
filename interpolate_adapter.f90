module interpolate
    use iso_fortran_env, only: dp => real64
    use bspline_module, only: bspline_1d, bspline_2d, bspline_3d

    implicit none

    type :: SplineData1D
        type(bspline_1d) :: s
    end type SplineData1D

    type :: SplineData2D
        type(bspline_2d) :: s
    end type SplineData2D

    type :: SplineData3D
        type(bspline_3d) :: s
    end type SplineData3D

contains

subroutine construct_splines_3d(x_min, x_max, fcn, order, periodic, spl)
    real(dp), intent(in) :: x_min(:), x_max(:), fcn(:,:,:)
    integer, intent(in) :: order(3)
    logical, intent(in) :: periodic(3)

    type(SplineData3D), intent(inout) :: spl

    integer :: kx, ky, kz, iflag, i, nx, ny, nz
    real(8) :: x(size(fcn, 1)), y(size(fcn, 2)), z(size(fcn, 3)), h

    nx = size(fcn, 1)
    ny = size(fcn, 2)
    nz = size(fcn, 3)

    kx = order(1) + 1
    ky = order(2) + 1
    kz = order(3) + 1

    h = (x_max(1) - x_min(1))/(nx-1)
    do i = 1, nx
        x(i) = x_min(1) + (i-1)*h
    end do

    h = (x_max(2) - x_min(2))/(ny-1)
    do i = 1, ny
        y(i) = x_min(2) + (i-1)*h
    end do

    h = (x_max(3) - x_min(3))/(nz-2)
    do i = 1, nz
        z(i) = x_min(3) + (i-1)*h
    end do

    call spl%s%initialize(x,y,z,fcn,kx,ky,kz,iflag,extrap=.True.)
end subroutine construct_splines_3d


subroutine evaluate_splines_3d(spl, x, y)
    type(SplineData3D), intent(inout) :: spl
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: y

    integer :: iflag

    call spl%s%evaluate(x(1),x(2),x(3),0,0,0,y,iflag)
end subroutine evaluate_splines_3d


subroutine evaluate_splines_3d_der2(spl, x, y, dy, d2y)
    type(SplineData3D), intent(inout) :: spl
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: y, dy(3), d2y(6)

    integer :: iflag

    call spl%s%evaluate(x(1),x(2),x(3),0,0,0,y,iflag)
    call spl%s%evaluate(x(1),x(2),x(3),1,0,0,dy(1),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),0,1,0,dy(2),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),0,0,1,dy(3),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),2,0,0,d2y(1),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),1,1,0,d2y(2),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),1,0,1,d2y(3),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),0,2,0,d2y(4),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),0,1,1,d2y(5),iflag)
    call spl%s%evaluate(x(1),x(2),x(3),0,0,2,d2y(6),iflag)
end subroutine evaluate_splines_3d_der2


subroutine destroy_splines_3d(spl)
    type(SplineData3D), intent(inout) :: spl
    call spl%s%destroy()
end subroutine destroy_splines_3d

end module interpolate
