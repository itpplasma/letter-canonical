# %%
import numpy as np
import matplotlib.pyplot as plt

shape = (64, 75, 50)

def doplot(data, title):
    # Plot 2D slice as contour
    plt.figure()
    plt.contour(data[:, 25, :])
    plt.colorbar()
    plt.title(title)

    plt.figure()
    plt.contour(data[:, :, 25])
    plt.colorbar()
    plt.title(title)

    plt.figure()
    plt.plot(data[:, 25, 25], 'o-')
    plt.title(title)


data = np.loadtxt("lam_phi.out").reshape(shape)
doplot(data, "lam_phi")

data = np.loadtxt("lam_phi_splined.out").reshape(shape)
doplot(data, "lam_phi_splined")

data = np.loadtxt("chi_gauge.out").reshape(shape)
doplot(data, "chi_gauge")

data = np.loadtxt("chi_gauge_splined.out").reshape(shape)
doplot(data, "chi_gauge_splined")

data = np.loadtxt("test_function.out").reshape(shape)
doplot(data, "test_function")

data = np.loadtxt("test_function_splined.out").reshape(shape)
doplot(data, "test_function_splined")


# %%
