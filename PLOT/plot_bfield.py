#%%
import numpy as np
import matplotlib.pyplot as plt

rmin = 75
rmax = 264.42281879194627
zmin = -150
zmax = 147.38193979933115

plt.figure()
data = np.loadtxt("fort.100")
for i in range(data.shape[0] - 1):
    if(np.sign(np.mod(data[i, 1], 2*np.pi) - np.pi) !=
       np.sign(np.mod(data[i+1, 1], 2*np.pi) - np.pi)):
        r = 0.5*(data[i, 0] + data[i+1, 0])
        z = 0.5*(data[i, 2] + data[i+1, 2])
        plt.plot(r, z, 'b,')


data = np.loadtxt("fort.101")
for i in range(data.shape[0] - 1):
    if(np.sign(np.mod(data[i, 1], 2*np.pi) - np.pi) !=
       np.sign(np.mod(data[i+1, 1], 2*np.pi) - np.pi)):
        r = 0.5*(data[i, 0] + data[i+1, 0])
        z = 0.5*(data[i, 2] + data[i+1, 2])
        plt.plot(r, z, 'r,')

plt.xlabel('R [cm]')
plt.ylabel('Z [cm]')
plt.xlim([rmin, rmax])
plt.ylim([zmin, zmax])

#%% 3D plot
from mpl_toolkits.mplot3d import Axes3D

tmax = 1000

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
data = np.loadtxt("fort.100")[:tmax,:]
X = data[:,0]*np.cos(data[:,1])
Y = data[:,0]*np.sin(data[:,1])
Z = data[:,2]
ax.plot(X, Y, Z)
data = np.loadtxt("fort.101")[:tmax,:]
X = data[:,0]*np.cos(data[:,1])
Y = data[:,0]*np.sin(data[:,1])
Z = data[:,2]
ax.plot(X, Y, Z, '--')
ax.set_xlim([-rmax, rmax])
ax.set_ylim([-rmax, rmax])
ax.set_zlim([zmin, zmax])

plt.show()

# %%
