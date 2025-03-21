#%%
%matplotlib tk
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
#%%

# Load data
data = np.loadtxt('../rk45_can.out')
R = data[:,0]
phi = data[:,1]
Z = data[:,2]
p = data[:,3]
vpar = data[:,4]

plt.plot(R, Z, '-')
#plt.xlim(1.35e7, 1.59e7)
#plt.ylim(-110, 100)

data = np.loadtxt('../rk45_cyl.out')
Rcyl = data[:,0]
phicyl = data[:,1]
Zcyl = data[:,2]
vpar_cyl = data[:,4]
plt.plot(Rcyl, Zcyl, '--')
# %%
plt.figure()
plt.plot(vpar, '-')
plt.plot(vpar_cyl, '--')
plt.xlabel("Timestep")
plt.ylabel("Parallel velocity")

# %%
plt.figure()
plt.plot(R, '-')
plt.plot(Rcyl, '--')
plt.xlabel("Timestep")
plt.ylabel("R(t)")

# %%
plt.figure()
plt.plot(phi)
plt.plot(phicyl)
plt.xlabel("Timestep")
plt.ylabel("phi(t)")

# %%
plt.figure()
plt.plot(Z, '-')
plt.plot(Zcyl, '--')
plt.xlabel("Timestep")
plt.ylabel("Z(t)")
# %%
