#%%
import numpy as np
import matplotlib.pyplot as plt

rmin = 75
rmax = 264.42281879194627
zmin = -150
zmax = 147.38193979933115

def plot_poincare_cuts(filename, linestyle):
    data = np.loadtxt(filename)
    for i in range(data.shape[0] - 1):
        if(np.sign(np.mod(data[i, 1], 2*np.pi) - np.pi) !=
        np.sign(np.mod(data[i+1, 1], 2*np.pi) - np.pi)):
            r = data[i, 0]
            z = data[i, 2]
            r = max(min([r, rmax]), rmin)
            z = max(min([z, zmax]), zmin)
            plt.plot(r, z, linestyle)
    plt.xlabel('R [cm]')
    plt.ylabel('Z [cm]')
    plt.axis('equal')


plt.figure()
plot_poincare_cuts("../fort.100", "b,")
plot_poincare_cuts("../fort.101", "r,")
plot_poincare_cuts("../fort.103", "g,")

#%%


#%% RZ plane plot
plt.figure()
data = np.loadtxt("../fort.100")
plt.plot(data[:, 0], data[:, 2], 'b,')
data = np.loadtxt("../fort.101")
plt.plot(data[:, 0], data[:, 2], 'r,')
plt.xlabel('R [cm]')
plt.ylabel('Z [cm]')
plt.axis('equal')


#%% 3D plot
from mpl_toolkits.mplot3d import Axes3D

tmax = 1500

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
data = np.loadtxt("../fort.100")[:tmax,:]
X = data[:,0]*np.cos(data[:,1])
Y = data[:,0]*np.sin(data[:,1])
Z = data[:,2]
ax.plot(X, Y, Z)
data = np.loadtxt("../fort.101")[:tmax,:]
X = data[:,0]*np.cos(data[:,1])
Y = data[:,0]*np.sin(data[:,1])
Z = data[:,2]
ax.plot(X, Y, Z, '--')
ax.set_xlim([-rmax, rmax])
ax.set_ylim([-rmax, rmax])
ax.set_zlim([zmin, zmax])

plt.show()

# %%
