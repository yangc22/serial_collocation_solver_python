import ex4test
import collocation_node
import numpy as np

'''
a = collocation_node.collocation_node(4, 1, 4, 3)
print(a.A)
a.R = np.eye(3)
print(a.R)
b = a.R + 1
print(b)

'''
bvp_dae = ex4test.bvp_dae()
print(bvp_dae.Y0.shape)
print(bvp_dae.Y0)
print(bvp_dae.Z0.shape)
print(bvp_dae.Z0)
print(bvp_dae.P0.shape)
print(bvp_dae.P0)