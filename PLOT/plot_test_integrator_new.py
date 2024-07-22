#%%
import numpy as np
import matplotlib.pyplot as plt

shape = (75, 64, 100)  # (Z, phi, R)

zmin = -150
zmax = 147
psimin = 12710026
psimax = 16377747

RR, ZZ = np.meshgrid(np.linspace(psimin, psimax, shape[2]),
    np.linspace(zmin, zmax, shape[0]))

#RZ plane plot of field
plt.figure()

data = np.loadtxt("../Aph_of_xc.out").reshape(shape)
plt.contour(RR, ZZ, data[:, 1, :], levels=20)
plt.xlabel("psi")
plt.ylabel("Z")

# Load data
data = np.loadtxt('../orbit.out')
R = data[:,0]
Z = data[:,1]
phi = data[:,2]
H = data[:,4]

plt.plot(R, Z, '-')
plt.xlim(psimin, psimax)
plt.ylim(zmin, zmax)

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
