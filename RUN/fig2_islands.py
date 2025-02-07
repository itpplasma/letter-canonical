#%%
import numpy as np
import matplotlib.pyplot as plt
from exportfig import exportpng

def load_banana(prefix):
    data = np.loadtxt(prefix + '_phi0minus.out')
    tau = data[:,0]
    z = data[:,1:]
    return tau, z

tau_midpoint, z_midpoint = load_banana('fig2_midpoint')    # 98021631 evals
tau_rk45_cyl, z_rk45_cyl = load_banana('fig2_rk45_cyl')    # 105489366 evals
tau_dop853_cyl, z_dop853_cyl = load_banana('fig2_dop853_cyl')  # 108608574 evals

plt.figure(figsize=(2.4,3.2))
plt.plot(z_dop853_cyl[:,0], z_dop853_cyl[:,2], ',', color='darkgray')
plt.plot(z_rk45_cyl[:,0], z_rk45_cyl[:,2], ',', color='lightgray')
plt.plot(z_midpoint[:,0], z_midpoint[:,2], ',', color='black')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')

plt.xlim([110, 240])
plt.ylim([-60, 70])

exportpng('fig2_islands')

# %%
