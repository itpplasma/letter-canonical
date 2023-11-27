"""
Christopher Albert, 2023

Does some checks and plots of the equilibrium + perturbation field interface
"""

# %%
import numpy as np

equilibrium_file = "../GORILLA/MHD_EQUILIBRIA/g_file_for_test"
convex_wall_file = "../GORILLA/MHD_EQUILIBRIA/convex_wall_for_test"

# %%
from my_little_magfie import my_little_magfie

my_little_magfie.init(equilibrium_file, convex_wall_file)

# %%

x = np.zeros(3)
B = np.zeros(3)

x[0] = 0.5
x[1] = 0.3
x[2] = 0.1

print(B)

my_little_magfie.eval_field_b(x, B)

print(B)
