from pathlib import Path
import math
import os
import subprocess
import sys
import tempfile
def write_inputs(root):
    n=33
    with (root/'synthetic.g').open('w') as f:
        f.write(f'{"synthetic circular equilibrium":48s}{0:4d}{n:4d}{n:4d}\n')
        def values(v):
            for i in range(0,len(v),5): f.write(''.join(f'{x:16.9E}' for x in v[i:i+5])+'\n')
        for row in [[2,2,2,1,0],[2,0,0,.3,2],[1e6,0,0,2,0],[0,0,.3,0,0]]: values(row)
        for row in [[4]*n,[0]*n,[0]*n,[0]*n]: values(row)
        values([.5*((1+2*i/(n-1)-2)**2+(-1+2*j/(n-1))**2) for j in range(n) for i in range(n)])
        values([1]*n)
        f.write('0 0\n')
    with (root/'wall.dat').open('w') as f:
        for i in range(64):
            a=2*math.pi*i/64
            f.write(f'{200+90*math.cos(a)} {90*math.sin(a)}\n')
    (root/'field.inp').write_text("0\n1\n1.0\n1\n0.99\n4\n'synthetic.g'\n''\n'wall.dat'\n''\n0\n0\n1\n")
    (root/'benchmark_gorilla.in').write_text('''&config
    magfie_type='tok', integrator_type='rk45', spatial_coordinates='cyl',
    velocity_coordinate='vpar', input_file_tok='field.inp',
    R0=210.0, phi0=0.0, Z0=12.0, vpar0=0.7,
    rtol=1d-9, dtau=100.0, ntau=250, energy_ev=3500.0,
    n_particles=8, poincare_shift=.false.
    /\n''')


def run(executable, root, threads):
    root.mkdir()
    write_inputs(root)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OMP_THREAD_LIMIT=str(threads))
    with (root / "run.log").open("w") as log:
        try:
            subprocess.run([executable, "benchmark_gorilla.in"], cwd=root, env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True, timeout=90)
        except (subprocess.SubprocessError, OSError):
            print((root / "run.log").read_text())
            raise
    results = {}
    for name in ("phi0plus", "phi0minus", "vpar0plus", "vpar0minus", "zplane"):
        rows = []
        for line in (root / ("benchmark_gorilla_" + name + ".out")).read_text().splitlines():
            parts = line.split()
            assert len(parts) == 7, "each callback record must include its orbit ID"
            orbit = int(parts[-1])
            values = tuple(float(value) for value in parts[:-1])
            assert 1 <= orbit <= 8 and all(math.isfinite(value) for value in values)
            rows.append((orbit, *values))
        results[name] = sorted(rows)
    assert {row[0] for row in results["phi0plus"]} == set(range(1, 9)), \
        "the synthetic run must exercise all eight orbits"
    return results


executable = str(Path(sys.argv[1]).resolve())
with tempfile.TemporaryDirectory(prefix="parallel-orbits-", dir=Path.cwd()) as scratch:
    root = Path(scratch)
    serial = run(executable, root / "serial", 1)
    for attempt in range(2):
        parallel = run(executable, root / ("parallel-" + str(attempt)), 4)
        assert parallel == serial, "per-orbit callback records differ between one and four threads"
print("Eight synthetic orbits agree exactly between one and four threads")
