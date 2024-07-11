#%%
import numpy as np
import matplotlib.pyplot as plt

# RZ plane plot of field
plt.figure()
data = np.loadtxt("../fort.6603")

R = data[:, 0]
Z = data[:, 1]
phi = data[:, 2]
pphi = data[:, 3]
f1 = data[:, -2]
f2 = data[:, -1]

plt.figure()
plt.plot(R, f1)

plt.figure()
plt.plot(R, f2)
# %%
