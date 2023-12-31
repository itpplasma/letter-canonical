# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fort.10")

data = data.reshape((75, 64, 50))

# Plot 2D slice as contour
plt.figure()
plt.contour(data[1, :, :])
plt.colorbar()


plt.figure()
plt.contour(data[:, :, 10])
plt.colorbar()

plt.figure()
plt.plot(data[10, :, 32], 'o-')

# %%
