import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fort.100"); plt.plot(data[:,0], data[:,1])
data = np.loadtxt("fort.101"); plt.plot(data[:,0], data[:,1])
