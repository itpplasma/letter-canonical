# %%
import numpy as np
import matplotlib.pyplot as plt

shape = (64, 75, 50)

# %% Plot delta_phi

data = np.loadtxt("lam_phi.out").reshape(shape)

# Plot 2D slice as contour
plt.figure()
plt.contour(data[25, :, :])
plt.colorbar()
plt.title("lam_phi")

plt.figure()
plt.contour(data[:, :, 25])
plt.colorbar()
plt.title("lam_phi")

plt.figure()
plt.plot(data[25, :, 25], 'o-')
plt.title("lam_phi")

# %% Plot gauge function chi

data = np.loadtxt("chi_gauge.out").reshape(shape)

# Plot 2D slice as contour
plt.figure()
plt.contour(data[25, :, :])
plt.colorbar()
plt.title("chi_gauge")

plt.figure()
plt.contour(data[:, :, 25])
plt.colorbar()
plt.title("chi_gauge")

plt.figure()
plt.plot(data[25, :, 25], 'o-')
plt.title("chi_gauge")

# %%
