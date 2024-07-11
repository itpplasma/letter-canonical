#%%
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('../euler1.out')

# Plot
R = data[:,0]
Z = data[:,1]
phi = data[:,2]
H = data[:,4]

plt.plot(R, Z, '-', label='Euler')
plt.xlim(130, 190)
plt.ylim(-30, 30)

# %%
plt.figure()
plt.plot(H, '-', label='Euler')
plt.xlim(0,1000)
plt.ylim(0,1000)
# %%
