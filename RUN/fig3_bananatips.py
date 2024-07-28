#%%
import numpy as np
import matplotlib.pyplot as plt

def load_banana(prefix):
    data = np.loadtxt(prefix + '_vpar0minus.out')
    tau = data[:,0]
    z = data[:,1:]
    return tau, z

tau_euler, z_euler = load_banana('fig3_midpoint')
tau_rk45_cyl, z_rk45_cyl = load_banana('fig3_rk45_cyl')
tau_rk45_can, z_rk45_can = load_banana('fig3_rk45_can')

plt.figure(figsize=(3.5,5))
plt.plot(z_rk45_cyl[:,0], z_rk45_cyl[:,2], ',', color='darkgray')
plt.plot(z_rk45_can[:,0], z_rk45_can[:,2], ',', color='lightgray')
plt.plot(z_euler[:,0], z_euler[:,2], ',', color='black')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')

#plt.xlim([125, 210])
#plt.ylim([-55, 65])

# %%
