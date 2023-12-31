"""
Christopher Albert, 2023

Does some checks and plots of the equilibrium + perturbation field interface
"""

# %%
import numpy as np
import matplotlib.pyplot as plt

equilibrium_file = "/proj/plasma/DATA/AUG/EQDSK/g_file_for_test"
perturbation_file = "/proj/plasma/DATA/AUG/COILS/MESH3D_31135odd/field_small.dat"
convex_wall_file = "convexwall_aug_test.dat"

convex_wall = np.loadtxt(convex_wall_file)
Rwall = convex_wall[:, 0]
Zwall = convex_wall[:, 1]
plt.plot(Rwall, Zwall, "k")

plt.axis("equal")

# %%
from my_little_magfie import my_little_magfie

my_little_magfie.init()

# %%

x = np.zeros(3)
B = np.zeros(3)
A = np.zeros(3)

x[0] = 100.0
x[1] = 0.0
x[2] = 0.0

print(B)

my_little_magfie.eval_field_b(x, B)

print(B)

# %%
nr = 75
nph = 1
nz = 150

Rmin = 75.0
Rmax = 267.0
Zmin = -150.0
Zmax = 150.0

R = np.linspace(Rmin, Rmax, nr)
P = np.linspace(0.0, 2 * np.pi, nph)
Z = np.linspace(Zmin, Zmax, nz)

RR, PP, ZZ = np.meshgrid(R, P, Z, indexing="ij")

BR = np.empty((nr, nph, nz))
BP = np.empty((nr, nph, nz))
BZ = np.empty((nr, nph, nz))

AR = np.empty((nr, nph, nz))
AP = np.empty((nr, nph, nz))
AZ = np.empty((nr, nph, nz))

for kr in range(nr):
    for kph in range(nph):
        for kz in range(nz):
            x[0] = RR[kr, kph, kz]
            x[1] = PP[kr, kph, kz]
            x[2] = ZZ[kr, kph, kz]
            my_little_magfie.eval_field_b_and_a(x, B, A)

            BR[kr, kph, kz] = B[0]
            BP[kr, kph, kz] = B[1]
            BZ[kr, kph, kz] = B[2]

            AR[kr, kph, kz] = A[0]
            AP[kr, kph, kz] = A[1]
            AZ[kr, kph, kz] = A[2]

# levels = np.linspace(-100.0, 100.0, 51)
levels = np.linspace(-30000.0, 10000.0, 51)

# %%

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


# %%

plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], AR[:, 0, :], levels=np.linspace(-2e3, 2e3, 51))
plt.plot(Rwall, Zwall, "k")
plt.title("AR")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/AR.pdf")


plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], AP[:, 0, :], levels=np.linspace(0, 3e5, 21))
plt.plot(Rwall, Zwall, "k")
plt.title("AP")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/AP.pdf")


plt.figure(figsize=(5, 6.4))
plt.contour(RR[:, 0, :], ZZ[:, 0, :], AZ[:, 0, :], levels=np.linspace(1e7, 1.7e7, 21))
plt.plot(Rwall, Zwall, "k")
plt.title("AZ")
plt.axis("equal")
plt.colorbar()
plt.savefig("PLOT/AZ.pdf")

# %%
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

R = np.linspace(Rmin, Rmax, nr)
Z = np.linspace(Zmin, Zmax, nz)

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


# %% Plot over phi
nph = 128

P = np.linspace(0.0, 2 * np.pi, nph)

x[0] = 0.5*(Rmin + Rmax)
x[2] = 0.5*(Zmin + Zmax)

BR = np.empty(nph)
BP = np.empty(nph)
BZ = np.empty(nph)

for kph in range(nph):
    x[1] = P[kph]
    my_little_magfie.eval_field_b(x, B)

    BR[kph] = B[0]
    BP[kph] = B[1]
    BZ[kph] = B[2]

Bmod = np.sqrt(BR**2 + BP**2 + BZ**2)

plt.figure()
plt.plot(P, BR)
plt.xlabel("phi")
plt.ylabel("BR")

plt.figure()
plt.plot(P, BP)
plt.xlabel("phi")
plt.ylabel("BP")

plt.figure()
plt.plot(P, BZ)
plt.xlabel("phi")
plt.ylabel("BZ")

plt.figure()
plt.plot(P, Bmod)
plt.xlabel("phi")
plt.ylabel("|B|")

# %% Load perturbation field directly for comparison

from pert_field import read_field_dat, get_grid

Br, Bp, Bz, bounds = read_field_dat(perturbation_file)

R, P, Z = get_grid(bounds, Br.shape)

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

# %% Do Poincare plot for field lines
from scipy.integrate import solve_ivp

num_cut = 256
steps_per_cut = 1
num_orbit = 4

Rstart = np.linspace(190, 206, num_orbit)

def fieldline_direction(t, x):
    B = np.zeros(3)

    my_little_magfie.eval_field_b(x, B)

    # Normalize by contravariant Bphi component for field line pitch
    B[0] = B[0] / (B[1] / x[0])
    B[2] = B[2] / (B[1] / x[0])
    B[1] = 1.0

    return B

plt.figure(figsize=(10, 12.8))
for k in range(num_orbit):
    xstart = np.array([Rstart[k], 0.0, 0.0])

    t = np.linspace(0.0, 2*np.pi*num_cut, steps_per_cut*num_cut + 1)
    sol = solve_ivp(fieldline_direction, [0.0, 2*np.pi*num_cut], xstart, t_eval=t,
                    method="RK45", rtol=1e-8, atol=1e-8)

    plt.plot(sol.y[0, ::steps_per_cut], sol.y[2, ::steps_per_cut], ".")

plt.xlim(Rmin, Rmax)
plt.ylim(Zmin, Zmax)
plt.plot(Rwall, Zwall, "k")
#plt.axis("equal")

# %%

B = np.zeros(3)
A = np.zeros(3)
x = np.array([Rstart[0], 0.0, 0.0])
my_little_magfie.eval_field_b_and_a(x, B, A)

print(B)
print(A)

# %%
