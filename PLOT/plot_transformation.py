# %%
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator
from plot_util import *

#%%
data = np.loadtxt("../lam_phi.out").reshape(shape)
doplot(data, "lam_phi")
plt.figure()
plt.plot(phi, phi+data[:, -1, -1], 'o-')
#%%
data = np.loadtxt("../lam_spl.out").reshape(shape_spl)
doplot(data, "lam_phi_splined")
#%%
data = np.loadtxt("../dlamdr_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dr_splined")
#%%
data = np.loadtxt("../dlamdz_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dz_splined")
#%%
data = np.loadtxt("../dlamdp_spl.out").reshape(shape_spl)
doplot(data, "dlam_phi_dphi_splined")

#%%
data = np.loadtxt("../chi_gauge.out").reshape(shape)
doplot(data, "chi_gauge")
#%%
data = np.loadtxt("../chi_spl.out").reshape(shape_spl)
doplot(data, "chi_gauge_splined")
#%%
data = np.loadtxt("../dchidr_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dr_splined")
#%%
data = np.loadtxt("../dchidz_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dz_splined")
#%%
data = np.loadtxt("../dchidp_spl.out").reshape(shape_spl)
doplot(data, "dchi_gauge_dphi_splined")

#%%
data = np.loadtxt("../A2_spl.out").reshape(shape_spl)
doplot(data, "A2_spl")
#%%
data = np.loadtxt("../A3_spl.out").reshape(shape_spl)
doplot(data, "A3_spl")

# %%
data = np.loadtxt("../psi_of_x.out").reshape(shape)
doplot(data, "psi_of_x")

# %%
data = np.loadtxt("../R_of_xc.out").reshape(shape)
doplot(data, "R_of_xc", rname="psi")
# %%
data = np.loadtxt("../Aph_of_xc.out").reshape(shape)
doplot(data, "Aph_of_xc", rname="psi")
# %%
data = np.loadtxt("../hph_of_xc.out").reshape(shape)
doplot(data, "hph_of_xc", rname="psi")
# %%
data = np.loadtxt("../hth_of_xc.out").reshape(shape)
doplot(data, "hth_of_xc", rname="psi")
# %%
data = np.loadtxt("../Bmod_of_xc.out").reshape(shape)
doplot(data, "Bmod_of_xc", rname="psi")

plt.show()

# %%
