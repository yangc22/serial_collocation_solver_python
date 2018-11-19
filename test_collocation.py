from collocation_solver import *
import numpy as np

if __name__ == '__main__':
    n_y = 4
    n_z = 1
    n_p = 4
    N = 101
    y0 = np.ones((N, n_y), dtype = np.float64)
    '''
    y0[0][0] = -1
    y0[1][0] = -2
    y0[2][0] = -3
    y0[3][0] = -4
    y0[4][0] = -5
    '''
    z0 = np.ones((N, n_z), dtype = np.float64)
    p0 = np.ones((n_p), dtype = np.float64)
    ti = 0
    tf = 1
    tspan = np.linspace(ti, tf, N, dtype = np.float64)
    # print(y0, z0, p0, tspan)
    collocation_solver(y0, z0, p0, tspan, n_y, n_z, n_p)