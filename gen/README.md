Regenerate the stable finite-segment integral using an installed
[FortSym](https://github.com/lazy-fortran/fortsym) prefix:

```sh
gfortran -I "$FORTSYM_PREFIX/include" gen/segment_log.f90 \
  -L "$FORTSYM_PREFIX/lib" -Wl,-rpath,"$FORTSYM_PREFIX/lib" -lfortsym \
  -o /tmp/generate_segment_log
/tmp/generate_segment_log
git diff --exit-code -- SRC/generated/segment_log.inc
```

Run from the repository root. Normal builds use the committed include and do
not require FortSym. The formula integrates `1/sqrt(rho**2+s**2)` between the
endpoint projections, where `rho` is the perpendicular distance and `s` the
coordinate along the segment. The test uses the independent axial geometry and
retains the existing tolerance.
