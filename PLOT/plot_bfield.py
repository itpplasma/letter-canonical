import numpy as np
import matplotlib.pyplot as plt

R0 = 165.0

plt.figure()
data = np.loadtxt("fort.100"); plt.plot(data[:,0], data[:,1])
data = np.loadtxt("fort.101"); plt.plot(data[:,0], data[:,1], '--')

# 3D plot
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
data = np.loadtxt("fort.100")
X = R0 + data[:,0]*np.cos(data[:,2])
Y = R0 + data[:,0]*np.sin(data[:,2])
Z = data[:,1]
ax.plot(X, Y, Z)
data = np.loadtxt("fort.101")
X = R0 + data[:,0]*np.cos(data[:,2])
Y = R0 + data[:,0]*np.sin(data[:,2])
Z = data[:,1]
ax.plot(X, Y, Z, '--')
plt.show()
