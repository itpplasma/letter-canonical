# %%
import numpy as np
import matplotlib.pyplot as plt

shape = (75, 64, 50)

# %% Plot delta_phi

data = np.loadtxt("delta_phi.out").reshape(shape)

# Plot 2D slice as contour
plt.figure()
plt.contour(data[25, :, :])
plt.colorbar()
plt.title("delta_phi")

plt.figure()
plt.contour(data[:, :, 25])
plt.colorbar()
plt.title("delta_phi")

plt.figure()
plt.plot(data[25, :, 25], 'o-')
plt.title("delta_phi")

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
