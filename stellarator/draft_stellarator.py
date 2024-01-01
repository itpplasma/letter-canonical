# %%
import numpy as np
import matplotlib.pyplot as plt
from simple_cyl import my_little_magfie

my_little_magfie.test()
my_little_magfie.init()

# %%
B = my_little_magfie.eval_field_b([1.0, 0.0, 0.0])
A = my_little_magfie.eval_field_a([1.0, 0.0, 0.0])

print(f"B = {B}")
print(f"A = {A}")
# %%
Rmin = 0.436
Rmax = 2.436
Zmin = -1.0
Zmax = 1.0

nr = 64
nph = 4
nz = 64

R = np.linspace(Rmin, Rmax, nr)
PH = np.linspace(0, 2 * np.pi*(nph-1)/nph, nph)
Z = np.linspace(Zmin, Zmax, nz)

RR, PP, ZZ = np.meshgrid(R, PH, Z, indexing='ij')

BB = np.zeros((nr, nph, nz, 3))
AA = np.zeros((nr, nph, nz, 3))

for kr in range(nr):
    for kph in range(nph):
        for kz in range(nz):
            r = RR[kr, kph, kz]
            ph = PP[kr, kph, kz]
            z = ZZ[kr, kph, kz]
            BB[kr, kph, kz] = my_little_magfie.eval_field_b([r, ph, z])
            AA[kr, kph, kz] = my_little_magfie.eval_field_a([r, ph, z])

# %%
kph = 0

plt.figure()
plt.quiver(RR[:, kph, :], ZZ[:, kph, :], BB[:, kph, :, 0], BB[:, kph, :, 2])

plt.figure()
plt.contourf(RR[:, kph, :], ZZ[:, kph, :], BB[:, kph, :, 1],
             levels=np.linspace(-2.0, 2.0, 100))
plt.colorbar()

plt.figure()
plt.quiver(RR[:, kph, :], ZZ[:, kph, :], AA[:, kph, :, 0], AA[:, kph, :, 2])

plt.figure()
plt.contourf(RR[:, kph, :], ZZ[:, kph, :], AA[:, kph, :, 1])

# %% Do Poincare plot for field lines
from scipy.integrate import solve_ivp

num_cut = 128
steps_per_cut = 1
num_orbit = 2

Rstart = np.linspace(1.5, 2.0, num_orbit)


plt.figure(figsize=(10, 10))
for k in range(num_orbit):
    xstart = np.array([Rstart[k], 0.0, 0.0])

    t = np.linspace(0.0, 2*np.pi*num_cut, steps_per_cut*num_cut + 1)
    sol = solve_ivp(my_little_magfie.fieldline_direction,
                    [0.0, 2*np.pi*num_cut], xstart, t_eval=t,
                    method="RK45", rtol=1e-6, atol=1e-6)

    plt.plot(sol.y[0, ::steps_per_cut], sol.y[2, ::steps_per_cut], ",")

plt.xlim(Rmin, Rmax)
plt.ylim(Zmin, Zmax)
#plt.axis("equal")

# %%
