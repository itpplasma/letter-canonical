program generate_segment_log
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(/), asinh
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SNIPPET
    use fortsym_string, only: str, chars
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: along, length, rho, roots(1)
    type(kernel_spec_t) :: spec
    integer :: unit

    call arena%init()
    along = sym(arena, 'along')
    length = sym(arena, 'L')
    rho = sym(arena, 'rho')
    ! Integral of 1/sqrt(rho**2 + s**2) between the two endpoint projections.
    roots(1) = asinh(along/rho) + asinh((length - along)/rho)
    spec%name = str('segment_log')
    allocate(spec%args(3))
    spec%args(1) = str('along')
    spec%args(2) = str('L')
    spec%args(3) = str('rho')
    spec%mode = KERNEL_SNIPPET
    spec%generator = str('gen/segment_log.f90')
    spec%regenerate_command = str('see gen/README.md')
    allocate(spec%outputs(1))
    spec%outputs(1) = str('log_term')
    open(newunit=unit, file='SRC/generated/segment_log.inc', status='replace')
    write(unit, '(a)') chars(emit_kernel(roots, spec))
    close(unit)
end program generate_segment_log
