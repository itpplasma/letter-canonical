module cut_callback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use callback_base, only: callback_t
    implicit none

    integer, parameter :: poly_order = 5
    integer, parameter :: buffer_size = poly_order + 1
    integer, parameter :: index_previous = buffer_size/2 - 1
    integer, parameter :: index_current = buffer_size/2

    type, extends(callback_t) :: cut_callback_t
        procedure(distance_p), nopass, pointer :: distance => null()
        procedure(event_p), nopass, pointer :: event => null()
        integer(8) :: counter = 0
        real(dp), allocatable :: buffer(:,:)
        logical :: is_initialized = .false.
        contains
        procedure :: init
        procedure :: execute
        procedure :: is_cut
        procedure :: trigger_event
        procedure :: values_at_cut
    end type cut_callback_t

    interface
        function distance_p(t, z)
            import dp
            real(dp) :: distance_p
            real(dp), intent(in) :: t, z(:)
        end function distance_p
    end interface

    interface
        subroutine event_p(t, z)
            import dp
            real(dp), intent(in) :: t, z(:)
        end subroutine event_p
    end interface

    contains

    subroutine init(self, n)
        class(cut_callback_t), intent(inout) :: self
        integer, intent(in) :: n

        allocate(self%buffer(n+2, buffer_size))
        self%buffer = 0d0
        self%is_initialized = .true.
    end subroutine init


    subroutine execute(self, t, z)
        class(cut_callback_t), intent(inout) :: self
        real(dp), intent(in) :: t
        real(dp), intent(in) :: z(:)

        if (.not. self%is_initialized) call self%init(size(z))

        call shift(self%buffer)
        call append(self%buffer, t, z, self%distance(t, z))

        self%counter = self%counter + 1
        if (self%counter < buffer_size) return

        if (is_cut(self)) call self%trigger_event

    end subroutine execute


    function is_cut(self)
        logical :: is_cut
        class(cut_callback_t), intent(in) :: self

        associate (b => self%buffer)
            if (self%distance(b(1, index_previous), b(2:, index_previous)) < 0d0 .and. &
                self%distance(b(1, index_current), b(2:, index_current)) >= 0d0) then
                is_cut = .true.
                return
            end if
        end associate
        is_cut = .false.
    end function is_cut


    subroutine values_at_cut(self, t, z)
        class(cut_callback_t), intent(in) :: self
        real(dp), intent(out) :: t
        real(dp), intent(out) :: z(:)

        call interpolate_at_cut(self%buffer, t, z)
    end subroutine values_at_cut


    subroutine trigger_event(self)
        class(cut_callback_t), intent(inout) :: self
        real(dp) :: t, z(size(self%buffer, 1) - 2)

        call self%values_at_cut(t, z)
        call self%event(t, z)
    end subroutine trigger_event


    subroutine interpolate_at_cut(buffer, t, z)
        use plag_coeff_sub

        real(dp), intent(in) :: buffer(:,:)
        real(dp), intent(out) :: t
        real(dp), intent(out) :: z(:)

        real(dp) :: coef(0:0, buffer_size)

        call plag_coeff(buffer_size, 0, 0d0, buffer(size(buffer, 1), :), coef)

        t = dot_product(buffer(1,:), coef(0,:))
        z = matmul(buffer(2:size(buffer, 1)-1,:), coef(0,:))
    end subroutine interpolate_at_cut


    subroutine shift(buffer)
        real(dp), intent(inout) :: buffer(:,:)
        integer :: i

        do i = 1, buffer_size-1
            buffer(:,i) = buffer(:,i+1)
        end do
    end subroutine shift


    subroutine append(buffer, t, z, distance)
        real(dp), intent(inout) :: buffer(:,:)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: z(:)
        real(dp), intent(in) :: distance

        buffer(1, buffer_size) = t
        buffer(2:size(buffer, 1) - 1, buffer_size) = z
        buffer(size(buffer, 1), buffer_size) = distance
    end subroutine append

end module cut_callback
