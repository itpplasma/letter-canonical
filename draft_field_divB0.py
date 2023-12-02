"""
Christopher Albert, 2023

Does some checks and plots of the equilibrium + perturbation field interface
"""

# %%
import numpy as np
import matplotlib.pyplot as plt

equilibrium_file = "../GORILLA/MHD_EQUILIBRIA/g_file_for_test"
perturbation_file = "/proj/plasma/DATA/AUG/COILS/MESH3D_31135odd/field.dat"
#convex_wall_file = "../GORILLA/MHD_EQUILIBRIA/convex_wall_for_test.dat"
#convex_wall_file = "convexwall_aug_rmp.dat"
convex_wall_file = "convexwall_aug_test.dat"

convex_wall = np.loadtxt(convex_wall_file)
Rwall = convex_wall[:, 0]
Zwall = convex_wall[:, 1]
plt.plot(Rwall, Zwall, "k")

plt.axis("equal")

# %%
from my_little_magfie import my_little_magfie

# %%

x = np.zeros(3)
B = np.zeros(3)

x[0] = 100.0
x[1] = 0.0
x[2] = 0.0

print(B)

my_little_magfie.eval_field_b(x, B)

print(B)

# %%
nr = 150
nph = 1
nz = 300

R = np.linspace(75.0, 267.0, nr)
P = np.linspace(0.0, 2 * np.pi, nph)
Z = np.linspace(-150.0, 150.0, nz)

RR, PP, ZZ = np.meshgrid(R, P, Z, indexing="ij")

BR = np.empty((nr, nph, nz))
BP = np.empty((nr, nph, nz))
BZ = np.empty((nr, nph, nz))

for kr in range(nr):
    for kph in range(nph):
        for kz in range(nz):
            x[0] = RR[kr, kph, kz]
            x[1] = PP[kr, kph, kz]
            x[2] = ZZ[kr, kph, kz]
            my_little_magfie.eval_field_b(x, B)

            BR[kr, kph, kz] = B[0]
            BP[kr, kph, kz] = B[1]
            BZ[kr, kph, kz] = B[2]

levels = np.linspace(-100.0, 100.0, 51)

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], BR[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BR")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/BR.pdf")

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], BP[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BP")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/BP.pdf")

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], BZ[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BZ")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/BZ.pdf")


plt.figure(figsize=(8, 8))
Bmod = np.sqrt(BR ** 2 + BP ** 2 + BZ ** 2)

plt.contourf(
    RR[:, 0, :],
    ZZ[:, 0, :],
    Bmod[:, 0, :],
)
plt.plot(Rwall, Zwall, "k")

plt.title("|B|")
plt.xlabel("R")
plt.ylabel("Z")
plt.axis("equal")

plt.colorbar()

plt.tight_layout()
plt.savefig("PLOT/Bmod.pdf")
# %% Plot B as a vector field

skip = 4

fig, ax = plt.subplots(figsize=(4, 8))

ax.quiver(
    RR[::skip, 0, ::skip],
    ZZ[::skip, 0, ::skip],
    BR[::skip, 0, ::skip],
    BZ[::skip, 0, ::skip]
)

ax.plot(Rwall, Zwall, "k")

ax.axis("equal")
ax.set_xlabel("R")
ax.set_ylabel("Z")

plt.savefig("PLOT/vector.pdf")
# %% Plot streamlines

fig, ax = plt.subplots(figsize=(4, 8))

ax.streamplot(
    RR[:, 0, :].T,
    ZZ[:, 0, :].T,
    BR[:, 0, :].T,
    BZ[:, 0, :].T
)

ax.axis("equal")
ax.set_xlabel("R")
ax.set_ylabel("Z")

plt.savefig("PLOT/streamlines.pdf")
# %%

nr = 1024
nz = 8

x[1] = 0.5
x[2] = 40.0

R = np.linspace(75.0, 267.0, nr)
Z = np.linspace(-150.0, 150.0, nz)

RR, ZZ = np.meshgrid(R, Z, indexing="ij")

BR = np.empty((nr, nz))
BP = np.empty((nr, nz))
BZ = np.empty((nr, nz))

for kr in range(nr):
    for kz in range(nz):
        x[0] = RR[kr, kz]
        x[2] = ZZ[kr, kz]
        my_little_magfie.eval_field_b(x, B)

        BR[kr, kz] = B[0]
        BP[kr, kz] = B[1]
        BZ[kr, kz] = B[2]

plt.figure(figsize=(8, 8))

for kz in range(nz):
    plt.plot(R, BR[:, kz], label=f"Z={Z[kz]}")

# %% Load perturbation field directly for compatison

# Open the file for reading
with open(perturbation_file, 'r') as iunit:
    # Read header information
    nr, nph, nz, _ = map(int, next(iunit).split())
    rmin, rmax = map(float, next(iunit).split())
    pmin, pmax = map(float, next(iunit).split())
    zmin, zmax = map(float, next(iunit).split())

    # Initialize arrays
    Br = np.zeros((nr, nph, nz), dtype=float)
    Bp = np.zeros((nr, nph, nz), dtype=float)
    Bz = np.zeros((nr, nph, nz), dtype=float)

    # Read data into arrays
    for i in range(nr):
        for j in range(nph):
            for k in range(nz):
                Br[i, j, k], Bp[i, j, k], Bz[i, j, k] = map(float, next(iunit).split())

# %% Plot
import matplotlib.cm as cm

R = np.linspace(75.0, 267.0, nr)
P = np.linspace(0.0, 2 * np.pi, nph)
Z = np.linspace(-150.0, 150.0, nz)

RR, PP, ZZ = np.meshgrid(R, P, Z, indexing="ij")

levels = np.linspace(-100.0, 100.0, 51)

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], Br[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BR")
plt.axis("equal")
plt.colorbar()

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], Bp[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BP")
plt.axis("equal")
plt.colorbar()

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], Bz[:, 0, :], levels=levels)
plt.plot(Rwall, Zwall, "k")
plt.title("BZ")
plt.axis("equal")
plt.colorbar()

# %%
