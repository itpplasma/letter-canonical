#%%
import numpy as np
import matplotlib.pyplot as plt

# RZ plane plot of field
plt.figure()
data = np.loadtxt("../fort.100")
plt.plot(data[:, 0], data[:, 2], 'b,')
data = np.loadtxt("../fort.101")
plt.plot(data[:, 0], data[:, 2], 'r,')
plt.xlabel('R [cm]')
plt.ylabel('Z [cm]')
plt.axis('equal')

# Load data
data = np.loadtxt('../orbit.out')
R = data[:,0]
Z = data[:,1]
phi = data[:,2]
H = data[:,4]

plt.plot(R, Z, '-')
plt.xlim(75, 250)
plt.ylim(-60, 60)

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
