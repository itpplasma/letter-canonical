#%%
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('../orbit.out')

# Plot
R0 = 1.0
r = data[:,0]
theta = data[:,1]
phi = data[:,2]

R = R0 + r*np.cos(theta)
Z = r*np.sin(theta)
plt.plot(R, Z, ',', label='Euler')

# %%
