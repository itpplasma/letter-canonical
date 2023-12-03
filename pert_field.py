import numpy as np

def read_field_dat(filename):
    """
    Read version 4 perturbation file field.dat
    as used in field_divB0.f90
    """
    # Open the file for reading
    with open(filename, 'r') as f:
        # Read header information
        nr, nph, nz, _ = map(int, next(f).split())
        rmin, rmax = map(float, next(f).split())
        pmin, pmax = map(float, next(f).split())
        zmin, zmax = map(float, next(f).split())

        # Initialize arrays
        Br = np.zeros((nr, nph, nz), dtype=float)
        Bp = np.zeros((nr, nph, nz), dtype=float)
        Bz = np.zeros((nr, nph, nz), dtype=float)

        # Read data into arrays
        for i in range(nr):
            for j in range(nph):
                for k in range(nz):
                    Br[i, j, k], Bp[i, j, k], Bz[i, j, k] = map(float, next(f).split())

        bounds = ((rmin, rmax), (pmin, pmax), (zmin, zmax))

        return Br, Bp, Bz, bounds

def write_field_dat(filename, Br, Bp, Bz, bounds):
    """
    Write version 4 perturbation file field.dat
    as used in field_divB0.f90
    """
    # Open the file for writing
    with open(filename, 'w') as f:
        # Write header information
        nr, nph, nz = Br.shape
        f.write(f"{nr} {nph} {nz} 4\n")
        f.write(f"{bounds[0][0]} {bounds[0][1]}\n")
        f.write(f"{bounds[1][0]} {bounds[1][1]}\n")
        f.write(f"{bounds[2][0]} {bounds[2][1]}\n")

        # Write data
        for i in range(nr):
            for j in range(nph):
                for k in range(nz):
                    f.write(f"{Br[i, j, k]} {Bp[i, j, k]} {Bz[i, j, k]}\n")

def get_grid(bounds, shape):
    """
    Return set of points in cylindrical coordinates for a tensor-product grid
    """
    rmin, rmax = bounds[0]
    pmin, pmax = bounds[1]
    zmin, zmax = bounds[2]

    nr, nph, nz = shape

    R = np.linspace(rmin, rmax, nr)
    P = pmin + np.arange(nph) * (pmax - pmin) / nph  # Don't duplicate 2pi
    Z = np.linspace(zmin, zmax, nz)

    return R, P, Z
