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
data = np.loadtxt('../euler1.out')
R = data[:,0]
Z = data[:,1]
phi = data[:,2]
H = data[:,4]

plt.plot(R, Z, '-', label='Euler')
plt.xlim(75, 250)
plt.ylim(-60, 60)

# %%
plt.figure()
plt.plot(H, '-', label='Euler')
# %%
