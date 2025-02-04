#%%
import numpy as np
import matplotlib.pyplot as plt
from exportfig import exportpng

def load_banana(prefix):
    data = np.loadtxt(prefix + '_phi0minus.out')
    tau = data[:,0]
    z = data[:,1:]
    return tau, z

tau_midpoint, z_midpoint = load_banana('fig2_midpoint')    #
tau_rk45_cyl, z_rk45_cyl = load_banana('fig2_rk45_cyl')    # 75140725 evals
tau_rk45_can, z_rk45_can = load_banana('fig2_dop853_cyl')  # 67438507 evals

plt.figure(figsize=(2.4,3.2))
plt.plot(z_rk45_can[:,0], z_rk45_can[:,2], ',', color='lightgray')
plt.plot(z_rk45_cyl[:,0], z_rk45_cyl[:,2], ',', color='darkgray')
plt.plot(z_midpoint[:,0], z_midpoint[:,2], ',', color='black')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')

plt.xlim([130, 140])
plt.ylim([-23, 12])

exportpng('fig2_islands')

# %%
