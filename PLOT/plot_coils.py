import sys
import numpy as np
import matplotlib.pyplot as plt

coil_file = sys.argv[1]
data = np.loadtxt(coil_file, skiprows=1)
coils = {
    'x': data[:, 0],
    'y': data[:, 1],
    'z': data[:, 2],
    'current': data[:, 3]
}

# Plot in 3D colored by current
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(coils['x'], coils['y'], coils['z'],
    c=np.abs(coils['current']), cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Coils')
plt.axis('equal')
plt.show()
