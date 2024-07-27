#%%
import numpy as np
import matplotlib.pyplot as plt

def load_data(prefix):
    data = np.loadtxt(prefix + '_pphi_H.out')
    pphi = data[:,0]
    H = data[:,1]
    return pphi, H

pphi_euler, H_euler = load_data('fig1_expl_impl_euler')
pphi_rk45_cyl, H_rk45_cyl = load_data('fig1_rk45_cyl')
pphi_rk45_can, H_rk45_can = load_data('fig1_rk45_can')

plt.figure()
plt.plot(H_euler, '-')
plt.plot(H_rk45_cyl, '-')
plt.plot(H_rk45_can, '-')
plt.xlabel('time step')
plt.ylabel('H')
plt.ylim([0.9, 1.1])


plt.figure()
plt.plot(pphi_euler, '-')
plt.plot(pphi_rk45_cyl, '-')
plt.plot(pphi_rk45_can, '-')
plt.xlabel('time step')
plt.ylabel('pphi')
plt.ylim([2e2, 2.4e2])

# %%
