import math
import numpy as np
from numba import jit, cuda, float32
import collocation_solver

@cuda.jit
def fKernel3D(d_f, d_x, d_y, d_z):
    i,j,k = cuda.grid(3)
    m,n,p = d_f.shape   
    if i < N:
        d_f[i,j,k] = collocation_solver.collocation_update(y[i], yDot[j], z[k])

def fArray3D(size_y, m, N, rk, sol):
    TPBX, TPBY, TPBZ = 8, 8, 8
    y = sol.get_y()
    yDot = sol.get_yDot()
    z = sol.get_z()
    y = cuda.to_device(y)
    d_yDot = cuda.to_device(y)
    d_z = cuda.to_device(z)
    F = cuda.device_array(shape = [N * m * (size_y + size_z) + N * size_y + size_p], dtype = np.float64)
    gridDims = (m+TPBX-1)/TPBX, (n+TPBY-1)/TPBY, (n+TPBZ-1)/TPBZ
    blockDims = TPBX, TPBY, TPBZ
    fKernel3D[gridDims, blockDims](d_f, d_x, d_y, d_z)

    return d_f.copy_to_host()