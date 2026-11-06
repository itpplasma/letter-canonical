program generate_segment_log
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(/), asinh
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SNIPPET
    use fortsym_string, only: str, chars, words
    implicit none
    type(expr_t) :: along, length, rho, roots(1)
    type(kernel_spec_t) :: spec
    integer :: unit

    along = sym('along')
    length = sym('L')
    rho = sym('rho')
    ! Integral of 1/sqrt(rho**2 + s**2) between the two endpoint projections.
    roots(1) = asinh(along/rho) + asinh((length - along)/rho)
    spec%name = str('segment_log')
    spec%args = words('along L rho')
    spec%mode = KERNEL_SNIPPET
    spec%generator = str('gen/segment_log.f90')
    spec%regenerate_command = str('see gen/README.md')
    spec%outputs = [str('log_term')]
    open(newunit=unit, file='SRC/generated/segment_log.inc', status='replace')
    write(unit, '(a)') chars(emit_kernel(roots, spec))
    close(unit)
end program generate_segment_log
