Paper on Overleaf: <https://www.overleaf.com/7388537417rqszybvwmksm#ccca0a>

```
make paper
```

to clone.

Fig. 1: Banana orbit crossing the separatrix in axisymmetric geometry. Symplectic schemes conserve toroidal momentum exactly and energy within bounds. Settings were such that all integrators used about 2000 field evaluations per bounce period. In this example, explicit-implicit Euler could still increase the step length and therefore reduce these evaluations by factor 3.9 until it reaches its stability limit.

Fig. 2: Islands of passing orbits
Callback output files keep time and the five phase-space coordinates in their
first six columns and append the orbit index as the seventh column. Parallel
records may be interleaved; group by orbit index before plotting a trajectory.
Each orbit starts with fresh callback history.

Cylindrical field evaluations serialize access to libneo's shared backend
state. The OpenMP region provides safe concurrent orbit state; it does not
promise a speedup for this serialized field path. The `parallel_orbits` CTest
case compares all callback records from eight synthetic orbits with one and
four threads, without experimental input data.
