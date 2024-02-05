# %%
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

shape = (64, 75, 50)
shape_spl = (63, 74, 49)

def doplot(data, title):
    # Plot 2D slice as contour
    plt.figure()
    plt.contour(data[:, 1, :])
    plt.colorbar()
    plt.title(title)

    plt.figure()
    plt.contour(data[:, :, 1])
    plt.colorbar()
    plt.title(title)

    plt.figure()
    plt.plot(data[:, 1, 1], 'o-')
    plt.plot(data[:, 2, 2], 'o-')
    plt.title(title)


data = np.loadtxt("lam_phi.out").reshape(shape)
doplot(data, "lam_phi")

data = np.loadtxt("lam_spl.out").reshape(shape_spl)
doplot(data, "lam_phi_splined")

data = np.loadtxt("chi_gauge.out").reshape(shape)
doplot(data, "chi_gauge")

data = np.loadtxt("chi_spl.out").reshape(shape_spl)
doplot(data, "chi_gauge_splined")


# %%
