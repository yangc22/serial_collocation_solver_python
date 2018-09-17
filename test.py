from collocation_coefficients import *
from collocation_node import *
import numpy as np
from ex4 import *

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

if __name__ == '__main__':
    bvp_dae = bvp_dae()
    y0 = np.asarray([1, 2, 3, 4], dtype = np.float64)
    y1 = np.asarray([5, 6, 7, 8], dtype = np.float64)
    z = np.asarray([1, 2, 3, 4, 5], dtype = np.float64)
    p = np.asarray([1, 2, 3], dtype = np.float64)
    alpha = 1
    f = np.zeros((4), dtype = np.float64)
    g = np.zeros((5), dtype = np.float64)
    r = np.zeros((7), dtype = np.float64)
    df = np.zeros((4, 12), dtype = np.float64)
    dg = np.zeros((5, 12), dtype = np.float64)
    dr = np.zeros((7, 11), dtype = np.float64)
    bvp_dae._abvp_f(y0, z, p, f)
    bvp_dae._abvp_g(y0, z, p, alpha, g)
    bvp_dae._abvp_r(y0, y1, p, r)
    bvp_dae._abvp_Df(y0, z, p, df)
    bvp_dae._abvp_Dg(y0, z, p, alpha, dg)
    bvp_dae._abvp_Dr(y0, y1, p, dr)
    print (f)
    print (g)
    print (r)
    print (df)
    print (dg)
    print (dr)
