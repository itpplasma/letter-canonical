# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("delta_phi.out")

data = data.reshape((75, 64, 50))

# Plot 2D slice as contour
plt.figure()
plt.contour(data[25, :, :])
plt.colorbar()


plt.figure()
plt.contour(data[:, :, 25])
plt.colorbar()

plt.figure()
plt.plot(data[25, :, 25], 'o-')

# %%
