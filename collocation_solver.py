from collocation_coefficients import *
from collocation_node import *
import numpy as np
from bvp_problem import *

'''
    input:
          ny - number of ODE variables
          nz - number of DAE variables
          np - number of parameter variables
          tspan - time span of the problem, with N time nodes
          y0 - initial estimate of ODE variables, N x ny matrix
          z0 - initial estimate of DAE variables, N x nz matrix
          p - initial estimate of parameter variables, np vector
'''
def collocation_solver(y0, z0, p0, tspan, n_y, n_z, n_p):
    size_y = n_y
    size_z = n_z
    size_p = n_p
    N = len(tspan) # number of time nodes
    m = 5 # number of coloocation points

    tol = 1e-6
    max_iter = 500
    max_linsearch = 20
    nodes_min = 3
    alpha = 1 # coontinuation parameter
    beta = 0.8 # scale factor

    rk = lobatto(m)
    bvp_dae = bvp_problem()

    sol = form_initial_input(size_y, size_z, size_p, m, N, tspan, y0, z0, p0, alpha, rk, bvp_dae)
    for i in range(N):
        print(sol[i].y, sol[i].z)
    for i in range(N):
        print(sol[i].tspan)

    tspan0, q0 = struct_to_vec(size_y, size_z, size_p, m, N, sol)
    sol = vec_to_struct(size_y, size_z, size_p, m, N, rk, tspan0, q0)
    #print(q0)
    # print(tspan0)
    print(sol[N - 1].y, sol[N - 1].p)
    F = F_q(size_y, size_z, size_p, m, N, rk, tspan0, q0, alpha)

'''
    form the initial input of collocation points from y0, z0, p0
'''
def form_initial_input(size_y, size_z, size_p, m, N, tspan, y0, z0, p0, alpha, rk, bvp_dae):
    sol = []
    c = rk.c
    for i in range(N - 1):
        node_i = collocation_node(size_y, size_z, size_p, m)
        node_i.set_y(y0[i][0 : size_y])
        if size_z > 0:
            node_i.set_z(z0[i][0 : size_z])
        tspan_i = np.zeros((m), dtype = np.float64)
        delta_t = tspan[i + 1] - tspan[i]
        node_i.set_delta_t(delta_t)
        for j in range(m):
            y_tmp = (1 - c[j]) * y0[i][0 : size_y] + c[j] * y0[i + 1][0 : size_y]
            z_tmp = (1 - c[j]) * z0[i][0 : size_z] + c[j] * z0[i + 1][0 : size_z]
            tspan_i[j] = tspan[i] + c[j] * delta_t
            y_dot = np.zeros((size_y), dtype = np.float64)
            bvp_dae._abvp_f(y_tmp, z_tmp, p0, y_dot)
            node_i.set_y_dot(y_dot, j)
            node_i.set_z_tilda(z_tmp, j)
        node_i.set_tspan(tspan_i)
        sol.append(node_i)

    node_N = collocation_node(size_y, size_z, size_p, m)
    node_N.set_y(y0[N - 1][0 : size_y])
    if size_z > 0:
        node_N.set_z(z0[i][0 : size_z])
    node_N.set_p(p0)
    node_N.set_tspan(tspan)
    sol.append(node_N)
    return sol

'''
    generate the vector form of the system variables from the structure form
'''
def struct_to_vec(size_y, size_z, size_p, m, N, sol):
    q = np.zeros((N * size_y + (N - 1) * m * (size_y + size_z) + size_p), dtype = np.float64)
    tspan = np.zeros((N), dtype = np.float64)
    for i in range(N):
        tspan[i] = sol[N - 1].tspan[i]
    for i in range(N - 1):
        start_index_y = i * (size_y + m * (size_y + size_z))
        for j in range(size_y):
            q[start_index_y + j] = sol[i].y[j]
        for j in range(m):
            start_index_ydot = start_index_y + size_y + j * (size_y + size_z)
            for k in range(size_y):
                q[start_index_ydot + k] = sol[i].y_dot[k][j]
            start_index_z = start_index_ydot + size_y
            for k in range(size_z):
                q[start_index_z + k] = sol[i].z_tilda[k][j]
    start_index_N = (N - 1) * size_y + (N - 1) * m * (size_y + size_z)
    for j in range(size_y):
        q[start_index_N + j] = sol[N - 1].y[j]
    start_index_p = start_index_N + size_y
    for j in range(size_p):
        q[start_index_p + j] = sol[N - 1].p[j]
    return tspan, q

'''
    generate the structure form of the system from the vector form
'''
def vec_to_struct(size_y, size_z, size_p, m, N, rk, tspan, q):
    c = rk.c
    sol = []
    for i in range(N - 1):
        node_i = collocation_node(size_y, size_z, size_p, m)
        tspan_i = np.zeros((m), dtype = np.float64)
        delta_t_i = tspan[i + 1] - tspan[i]
        node_i.set_delta_t(delta_t_i)
        start_index_y = i * (size_y + m * (size_y + size_z))
        end_index_y = start_index_y + size_y
        node_i.set_y(q[start_index_y : end_index_y])
        for j in range(m):
            tspan_i[j] = tspan[i] + c[j] * delta_t_i
            start_index_y_collocation = end_index_y + j * (size_y + size_z)
            end_index_y_collocation = start_index_y_collocation + size_y
            start_index_z_collocation = end_index_y_collocation
            end_index_z_collocation = start_index_z_collocation + size_z
            node_i.set_y_dot(q[start_index_y_collocation : end_index_y_collocation], j)
            node_i.set_z_tilda(q[start_index_z_collocation : end_index_z_collocation], j)
        node_i.set_tspan(tspan_i)
        sol.append(node_i)
    node_N = collocation_node(size_y, size_z, size_p, m)
    node_N.set_tspan(tspan)
    start_index_y = (N - 1) * size_y + (N - 1) * m * (size_y + size_z)
    end_index_y = start_index_y + size_y
    start_index_p = end_index_y
    end_index_p = start_index_p + size_p
    node_N.set_y(q[start_index_y : end_index_y])
    node_N.set_p(q[start_index_p : end_index_p])
    sol.append(node_N)
    return sol

'''
    calculate the value of ODE variables on each collocation point
'''
def  collocation_update(size_y, m, N, rk, sol):
    a = rk.A
    for i in range(N - 1):
        delta_t_i = sol[i].delta_t
        for j in range(m):
            sum_j = np.zeros((size_y), dtype = np.float64)
            for k in range(m):
                for l in range(size_y):
                    sum_j[l] += a[j][k] * sol[i].y_dot[l][k]
            '''
            y_tilda_j = np.zeros((size_y), dtype = float64)
            for l in range(size_y):
                y_tilda_j[l] = sol[i].y[l] + delta_t * sum_j[l]
            '''
            y_tilda_j = sol[i].y + delta_t_i * sum_j
            sol[i].set_y_tilda(y_tilda_j, j)

'''
    calculate the residual of the system
'''
def F_q(size_y, size_z, size_p, m, N, rk, tspan, q, alpha):
    F = np.zeros((N * size_y + (N - 1) * m * (size_y + size_z) + size_p), dtype = np.float64)
    sol = vec_to_struct(size_y, size_z, size_p, m, N, rk, tspan, q)
    collocation_update(size_y, m, N, rk, sol)
    return F