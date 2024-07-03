# %%
import numpy as np
from plot_util import *

data = np.loadtxt("../Aphi_of_xc.out").reshape(shape)
doplot3d(data, "Aphi_of_xc", "𝜓")

# %%
