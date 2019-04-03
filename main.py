from collocation_coefficients import *
from collocation_node import *
import numpy as np
import bvp_problem
import math
import form_initial_inpu_serial

print('hello')
if __name__ == '__main__':
    print('hi')
    bvp_dae = bvp_problem.bvp_dae()

    size_y = bvp_dae.size_y
    size_z = bvp_dae.size_z
    size_p = bvp_dae.size_p
    size_inequality = bvp_dae.size_inequality
    N = 3  # number of time nodes
    tspan = np.linspace(0, 1, N)
    y0 = np.arange(N * size_y).reshape((N, size_y))
    z0 = np.arange(N * size_z).reshape((N, size_z))
    p0 = np.ones(size_p)
    m = 3  # number of coloocation points
    rk = lobatto(m)

    sol = form_initial_inpu_serial.form_initial_input(bvp_dae, size_y, size_z, size_p, m, N, tspan, y0, z0, p0, alpha, rk)