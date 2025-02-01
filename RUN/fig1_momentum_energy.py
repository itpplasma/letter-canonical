#%%
import numpy as np
import matplotlib.pyplot as plt
from exportfig import exportpng

tau_bounce = 8373.8  # From vpar0_plus output

def load_data(prefix):
    data = np.loadtxt(prefix + '_pphi_H.out')
    tau = data[:,0]
    pphi = data[:,1]
    H = data[:,2]
    return tau, pphi, H

tau_midpoint, pphi_midpoint, H_midpoint = load_data('fig1_midpoint')
tau_rk45_cyl, pphi_rk45_cyl, H_rk45_cyl = load_data('fig1_rk45_cyl')
tau_rk45_can, pphi_rk45_can, H_rk45_can = load_data('fig1_dop853_cyl')

plt.figure()
plt.plot(tau_rk45_cyl/tau_bounce, H_rk45_cyl/H_rk45_cyl[0], ',', color='darkgray')
plt.plot(tau_rk45_can/tau_bounce, H_rk45_can/H_rk45_can[0], ',', color='lightgray')
plt.plot(tau_midpoint/tau_bounce, H_midpoint/H_midpoint[0], ',', color='black')
plt.xlabel(r'number of bounce periods $t/\tau_{\mathrm{b}}$')
plt.ylabel(r'normalized energy $H$')
plt.ylim([0.9, 1.1])
exportpng('fig1b_energy')

plt.figure()
plt.plot(tau_rk45_cyl/tau_bounce, pphi_rk45_cyl/pphi_rk45_cyl[0], ',', color='darkgray')
plt.plot(tau_rk45_can/tau_bounce, pphi_rk45_can/pphi_rk45_can[0], ',', color='lightgray')
plt.plot(tau_midpoint/tau_bounce, pphi_midpoint/pphi_midpoint[0], ',', color='black')
plt.xlabel(r'number of bounce periods $t/\tau_{\mathrm{b}}$')
plt.ylabel(r'normalized toroidal momentum $p_\varphi$')
plt.ylim([0.9, 1.1])
exportpng('fig1c_pphi')

# %%

def load_banana(prefix):
    data = np.loadtxt(prefix + '.out')
    tau = data[:,0]
    z = data[:,1:]
    return tau, z

tau_midpoint, z_midpoint = load_banana('fig1_midpoint')
tau_rk45_cyl, z_rk45_cyl = load_banana('fig1_rk45_cyl')
tau_dop853_cyl, z_dop853_cyl = load_banana('fig1_dop853_cyl')

plt.figure(figsize=(2.4,3.2))
plt.plot(z_rk45_cyl[:,0], z_rk45_cyl[:,2], ',', color='darkgray')
plt.plot(z_dop853_cyl[:,0], z_dop853_cyl[:,2], ',', color='lightgray')
plt.plot(z_midpoint[:,0], z_midpoint[:,2], ',', color='black')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')

plt.xticks([100, 150, 200, 250])
plt.xlim([100, 250])
plt.ylim([-132, 100])

exportpng('fig1a_banana')
# %%
