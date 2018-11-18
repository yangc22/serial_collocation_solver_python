import bvp_problem
import numpy as np

bvp_dae = bvp_problem.bvp_dae()
size_y = 4
size_z = 1
size_p = 4
y = np.ones((size_y))
z = np.ones((size_z))
p = np.ones((size_p))

Dh = np.zeros((size_y, (size_y + size_z + size_p)), dtype = np.float64)
print(Dh)
bvp_dae._abvp_Df(y, z, p, Dh)
print(Dh)