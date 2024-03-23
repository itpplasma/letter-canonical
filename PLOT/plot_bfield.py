import numpy as np
import matplotlib.pyplot as plt

rmin = 75
rmax = 264.42281879194627
zmin = -150
zmax = 147.38193979933115

plt.figure()
data = np.loadtxt("fort.100"); plt.plot(data[:,0], data[:,1])
data = np.loadtxt("fort.101"); plt.plot(data[:,0], data[:,1], '--')
plt.xlabel('R [cm]')
plt.ylabel('Z [cm]')
plt.xlim([rmin, rmax])
plt.ylim([zmin, zmax])

# 3D plot
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
data = np.loadtxt("fort.100")
X = data[:,0]*np.cos(data[:,2])
Y = data[:,0]*np.sin(data[:,2])
Z = data[:,1]
ax.plot(X, Y, Z)
data = np.loadtxt("fort.101")
X = data[:,0]*np.cos(-data[:,2])
Y = data[:,0]*np.sin(-data[:,2])
Z = data[:,1]
ax.plot(X, Y, Z, '--')
ax.set_xlim([-rmax, rmax])
ax.set_ylim([-rmax, rmax])
ax.set_zlim([zmin, zmax])

plt.show()
