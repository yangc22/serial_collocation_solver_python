from collocation_coefficients import *
from collocation_node import *
import numpy as np

'''
n = collocation_node(4,4,1,3);
y = np.zeros((4))
z = np.ones((4))
pp = [1]
p = np.asarray(pp)
print(p)

n.set_y(y)
n.set_z(z)
n.set_p(p)
print(n.y, n.z, n.p)
y[0] = 10
z[0] = 10
p[0] = 10
print(n.y, n.z, n.p)

lo = lobatto(4)
print(lo.A)

nodes = []
n1 = collocation_node(4, 4, 1, 1)
n2 = collocation_node(4, 4, 1, 1)
nodes.append(n1)
nodes.append(n2)
print(nodes[1].y)

yy = [1, 2, 3, 4]
zz = [1]
p = [1]
alpha = [1]
h_ode = ocp.ODE_h(yy, zz, p, alpha)
print(h_ode)
'''

def helper(y):
    y[0] = -1

if __name__ == '__main__':
    y = np.zeros((4), dtype = np.float64)
    print(y)
    helper(y)
    print(y)