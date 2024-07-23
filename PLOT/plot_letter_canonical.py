#%%
import numpy as np
import matplotlib.pyplot as plt

shape = (75, 64, 100)  # (Z, phi, R)
shape_spl = (74, 63, 49)

zmin = -150.0
zmax = 147.38193979933115
rmin = 75.0
rmax = 264.42281879194627

RR, ZZ = np.meshgrid(np.linspace(rmin, rmax, shape_spl[2]),
                     np.linspace(zmin, zmax, shape_spl[0]))

#RZ plane plot of field
plt.figure()

data = np.loadtxt("../A2_spl.out").reshape(shape_spl)
plt.contour(RR, ZZ, data[:, 1, :], levels=19)
plt.xlabel("R")
plt.ylabel("Z")


# Load data
data = np.loadtxt('../letter_canonical.out')
R = data[:,0]
phi = data[:,1]
Z = data[:,2]
vpar = data[:,3]

plt.plot(R, Z, '-')
#plt.xlim(1.35e7, 1.59e7)
#plt.ylim(-110, 100)

# %%
plt.figure()
plt.plot(H[2:], '-')
plt.xlabel("Timestep")
plt.ylabel("Hamiltonian H")
plt.ylim(H[2:-1].min()*0.999, H[2:-1].max()*1.001)
# %%

plt.figure()
plt.plot(Z)
plt.xlabel("Timestep")
plt.ylabel("Z(t)")
plt.ylim(Z[2:-1].min()*1.001, Z[2:-1].max()*0.999)
# %%


plt.figure()
plt.plot(phi)
plt.xlabel("Timestep")
plt.ylabel("phi(t)")
plt.ylim(phi[2:-1].min()*1.001, phi[2:-1].max()*0.999)
# %%
