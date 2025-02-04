#%%
import numpy as np
import matplotlib.pyplot as plt
from exportfig import exportpng

def load_banana(prefix):
    data = np.loadtxt(prefix + '_vpar0minus.out')
    tau = data[:,0]
    z = data[:,1:]
    return tau, z

tau_euler, z_euler = load_banana('fig3_midpoint')          # 44965629 evaluations
tau_rk45_cyl, z_rk45_cyl = load_banana('fig3_rk45_cyl')    # 47806918 evaluations
tau_dop853_cyl, z_dop853_cyl = load_banana('fig3_dop853_cyl')  # 42900006 evaluations

plt.figure(figsize=(2.4,3.2))
plt.plot(z_rk45_cyl[:,0], z_rk45_cyl[:,2], ',', color='darkgray')
plt.plot(z_dop853_cyl[:,0], z_dop853_cyl[:,2], ',', color='lightgray')
plt.plot(z_euler[:,0], z_euler[:,2], ',', color='black')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')

plt.xlim([153.6, 154.2])
plt.ylim([53, 61.5])
exportpng('fig3_bananatips')

# %%
