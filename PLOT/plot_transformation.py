# %%
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

shape = (75, 64, 100)
shape_spl = (74, 63, 49)

phi = np.linspace(0, 2 * np.pi, shape[0])

def doplot(data, title):
    # Plot 2D slice as contour
    try:
        plt.figure()
        plt.contour(data[:, :, 1])
        plt.colorbar()
        plt.xlabel("r")
        plt.ylabel("phi")
        plt.title(title)
    except:
        pass

    try:
        plt.figure()
        plt.contour(data[:, 1, :])
        plt.colorbar()
        plt.xlabel("z")
        plt.ylabel("phi")
        plt.title(title)
    except:
        pass

    try:
        plt.figure()
        plt.plot(data[1, :, 1], 'o-')
        plt.plot(data[2, :, 2], 'o-')
        plt.xlabel("phi")
        plt.title(title)
    except:
        pass

#%%
data = np.loadtxt("lam_phi.out").reshape(shape)
doplot(data, "lam_phi")
plt.figure()
plt.plot(phi, phi+data[:, -1, -1], 'o-')
#%%
data = np.loadtxt("lam_spl.out").reshape(shape_spl)
doplot(data, "lam_phi_splined")
#%%
data = np.loadtxt("dlamdr_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dr_splined")
#%%
data = np.loadtxt("dlamdz_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dz_splined")
#%%
data = np.loadtxt("dlamdp_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dphi_splined")

#%%
data = np.loadtxt("chi_gauge.out").reshape(shape)
doplot(data, "chi_gauge")
#%%
data = np.loadtxt("chi_spl.out").reshape(shape_spl)
doplot(data, "chi_gauge_splined")
#%%
data = np.loadtxt("dchidr_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dr_splined")
#%%
data = np.loadtxt("dchidz_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dz_splined")
#%%
data = np.loadtxt("dchidp_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dphi_splined")

#%%
data = np.loadtxt("A1_spl.out").reshape(shape_spl)
doplot(data, "A1_spl")
#%%
data = np.loadtxt("A2_spl.out").reshape(shape_spl)
doplot(data, "A2_spl")
#%%
data = np.loadtxt("A3_spl.out").reshape(shape_spl)
doplot(data, "A3_spl")

# %%

plt.show()
