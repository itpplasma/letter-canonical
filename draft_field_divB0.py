"""
Christopher Albert, 2023

Does some checks and plots of the equilibrium + perturbation field interface
"""

# %%
import numpy as np

equilibrium_file = "../GORILLA/MHD_EQUILIBRIA/g_file_for_test"
perturbation_file = "../GORILLA/MHD_EQUILIBRIA/netcdf_file_for_test.nc"
convex_wall_file = "../GORILLA/MHD_EQUILIBRIA/convex_wall_for_test"

# %%
import field_gorilla
