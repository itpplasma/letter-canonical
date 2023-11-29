"""
Christopher Albert, 2023

Does some checks and plots of the equilibrium + perturbation field interface
"""

# %%
import numpy as np
import matplotlib.pyplot as plt

equilibrium_file = "../GORILLA/MHD_EQUILIBRIA/g_file_for_test"
convex_wall_file = "../GORILLA/MHD_EQUILIBRIA/convex_wall_for_test.dat"

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
nr = 64
nph = 1
nz = 128

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

# %% Plot

fig, ax = plt.subplots(1, 3, figsize=(12, 4))

ax[0].contourf(RR[:, 0, :], ZZ[:, 0, :], BR[:, 0, :])
ax[0].set_title("BR")
ax[1].contourf(RR[:, 0, :], ZZ[:, 0, :], BP[:, 0, :])
ax[1].set_title("BP")
ax[2].contourf(RR[:, 0, :], ZZ[:, 0, :], BZ[:, 0, :])
ax[2].set_title("BZ")

plt.savefig("PLOT/Bfield.pdf")

for k in range(3):
    ax[k].set_xlabel("R")
    ax[k].set_ylabel("Z")

fig.tight_layout()

plt.figure()
plt.contourf(
    RR[:, 0, :],
    ZZ[:, 0, :],
    np.sqrt(BR[:, 0, :] ** 2 + BP[:, 0, :] ** 2 + BZ[:, 0, :] ** 2),
)
plt.title("|B|")
plt.xlabel("R")
plt.ylabel("Z")

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
