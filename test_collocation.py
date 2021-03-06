from collocation_solver import *
from bvp_problem import bvp_dae
import numpy as np

if __name__ == '__main__':
    '''
    n_y = 4
    n_z = 1
    n_p = 4
    N = 101
    y0 = np.ones((N, n_y), dtype = np.float64)
    z0 = np.ones((N, n_z), dtype = np.float64)
    p0 = np.ones((n_p), dtype = np.float64)
    ti = 0
    tf = 1
    tspan = np.linspace(ti, tf, N, dtype = np.float64)
    # print(y0, z0, p0, tspan)
    collocation_solver(y0, z0, p0, tspan, n_y, n_z, n_p)

    '''
    n_y = 4
    n_z = 5
    n_p = 3
    N = 101
    y0 = np.ones((N, n_y), dtype=np.float64)
    z0 = np.ones((N, n_z), dtype=np.float64)
    p0 = np.ones(n_p, dtype=np.float64)
    ti = 0
    tf = 1
    tspan = np.linspace(ti, tf, N, dtype=np.float64)
    # print(y0, z0, p0, tspan)
    # collocation_solver(y0, z0, p0, tspan, n_y, n_z, n_p)
    bvpDae = bvp_dae()
    collocation_solver(bvpDae)
