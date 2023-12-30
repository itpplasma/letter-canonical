# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fort.10")

data = data.reshape((75, 32, 50))

# Plot 2D slice as contour
plt.figure()
plt.contour(data[:, :, 10])
plt.colorbar()

# %%
