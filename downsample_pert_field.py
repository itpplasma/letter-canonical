# %%
import numpy as np

from pert_field import read_field_dat, write_field_dat, get_grid

perturbation_file = "/proj/plasma/DATA/AUG/COILS/MESH3D_31135odd/field.dat"

Br, Bp, Bz, bounds = read_field_dat(perturbation_file)

print(f"Original grid size: {Br.shape}")
print(f"Original bounds: {bounds}")

R, P, Z = get_grid(bounds, Br.shape)

# %% Downsample
nr_new = 50
np_new = 32
nz_new = 75

nr_factor = Br.shape[0] // nr_new
np_factor = Br.shape[1] // np_new
nz_factor = Br.shape[2] // nz_new

Br_new = Br[::nr_factor, ::np_factor, ::nz_factor]
Bp_new = Bp[::nr_factor, ::np_factor, ::nz_factor]
Bz_new = Bz[::nr_factor, ::np_factor, ::nz_factor]

R_new = R[::nr_factor]
P_new = P[::np_factor]
Z_new = Z[::nz_factor]

bounds_new = (
    (R_new.min(), R_new.max()),
    (P_new.min(), 2 * np.pi),
    (Z_new.min(), Z_new.max()),
)
# %% Write
write_field_dat("field_test.dat", Br_new, Bp_new, Bz_new, bounds_new)

# %%
