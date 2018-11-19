from collocation_coefficients import *
from collocation_node import *
import numpy as np
import bvp_problem
import math
import matplotlib.pyplot as plt

'''
    original MATLAB code sol(N).f_b is changed to sol[N - 1].f_N
'''

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
    m = 3 # number of coloocation points

    tol = 1e-6
    max_iter = 500
    max_linsearch = 20
    nodes_min = 3
    alpha = 1 # coontinuation parameter
    alpha_final = 1
    beta = 0.8 # scale factor

    rk = lobatto(m)
    bvp_dae = bvp_problem.bvp_dae()

    sol = form_initial_input(size_y, size_z, size_p, m, N, tspan, y0, z0, p0, alpha, rk, bvp_dae)

    tspan0, q0 = struct_to_vec(size_y, size_z, size_p, m, N, sol)

    '''
    sol = vec_to_struct(size_y, size_z, size_p, m, N, rk, tspan0, q0)

    F, sol = F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan0, q0, alpha)

    # M = f_d_jacobian(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan0, q0, alpha)
    Jacobian_construct(bvp_dae, size_y, size_z, size_p, m, N, rk, alpha, sol)

    for i in range(N):
        print("i: ", i)
        print("J: ", sol[i].J)
        print("V: ", sol[i].V)
        print("D: ", sol[i].D)
        print("W: ", sol[i].W)
        print("B: ", sol[i].B)
        print("V_N: ", sol[i].V_N)
        print("A: ", sol[i].A)
        print("C: ", sol[i].C)
        print("H: ", sol[i].H)
        print("b: ", sol[i].b)
        print("b_N: ", sol[i].b_N)

    qr_decomposition(size_y,size_p, N, sol)
    print(sol[0].R)
    print(sol[0].E)
    print(sol[0].G)
    print(sol[0].K)
    print(sol[0].d)
    print(sol[N - 2].d)
    print(sol[N - 1].d)
    print(sol[N - 1].dp)
    print(sol[N - 1].Rp)

    backward_substitution(N, sol)
    for i in range(N):
        print(sol[i].delta_y)
        print(sol[i].delta_k)

    delta_q = get_delta_q(size_y, size_z, size_p, m, N, sol)
    print(delta_q)
    '''

    for alphacal in range(max_iter):
        for iter_time in range(max_iter):
            F_q0, sol = F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan0, q0, alpha)
            norm_F_q0 = np.linalg.norm(F_q0, np.inf)
            if (norm_F_q0 < tol):
                print("solution found")
                break
            Jacobian_construct(bvp_dae, size_y, size_z, size_p, m, N, rk, alpha, sol)
            # M = f_d_jacobian(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan0, q0, alpha)
            # delta_q0 = np.linalg.solve(M, -F_q0)
            qr_decomposition(size_y,size_p, N, sol)
            backward_substitution(N, sol)
            delta_q0 = get_delta_q(size_y, size_z, size_p, m, N, sol)
            norm_delta_q0 = np.linalg.norm(delta_q0, np.inf)
            
            if (norm_delta_q0 < tol):
                print("solution found")
                break
            print("start lin search")
            # print(delta_q0)
            alpha0 = 1
            for i in range(max_linsearch):
                q_new = q0 + alpha0 * delta_q0
                F_q_new, _ = F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan0, q_new, alpha)
                norm_F_q_new = np.linalg.norm(F_q_new, np.inf)
                if (norm_F_q_new < norm_F_q0):
                    q0 = q_new
                    break
                alpha0 /= 2
        print("final solution found")
        if (alpha <= alpha_final):
            break
    print("done")


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
def F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan, q, alpha):
    F = np.zeros((N * size_y + (N - 1) * m * (size_y + size_z) + size_p), dtype = np.float64)
    sol = vec_to_struct(size_y, size_z, size_p, m, N, rk, tspan, q)
    collocation_update(size_y, m, N, rk, sol)
    p = sol[N - 1].p
    b = rk.b
    y = np.zeros((size_y, N), dtype = np.float64)
    for i in range(N):
        for j in range(size_y):
            y[j][i] = sol[i].y[j]
    for i in range(N - 1):
        f_a = np.zeros((m * (size_y + size_z)), dtype = np.float64)
        for j in range(m):
            r_h = np.zeros((size_y), dtype = np.float64)
            r_g = np.zeros((size_z), dtype = np.float64)
            y_tmp = np.zeros((size_y), dtype = np.float64)
            z_tmp = np.zeros((size_z), dtype = np.float64)
            for k in range(size_y):
                y_tmp[k] = sol[i].y_tilda[k][j]
            for k in range(size_z):
                z_tmp[k] = sol[i].z_tilda[k][j]
            bvp_dae._abvp_f(y_tmp, z_tmp, p, r_h)
            bvp_dae._abvp_g(y_tmp, z_tmp, p, alpha, r_g)
            for k in range(size_y):
                r_h[k] -= sol[i].y_dot[k][j]
            start_index_y = j * (size_y + size_z)
            for k in range(size_y):
                f_a[start_index_y + k] = r_h[k]
            start_index_z = start_index_y + size_y
            for k in range(size_z):
                f_a[start_index_z + k] = r_g[k]
        sol[i].set_f_a(f_a)
        sum_i = np.zeros((size_y), dtype = np.float64)
        for j in range(m):
            for k in range(size_y):
                sum_i[k] += b[j] * sol[i].y_dot[k][j]
        y_next = np.zeros((size_y), dtype = np.float64)
        y_cur = np.zeros((size_y), dtype = np.float64)
        for j in range(size_y):
            y_next[j] = y[j][i + 1] 
            y_cur[j] = y[j][i]
        delta_t = sol[i].delta_t
        r_y = np.zeros((size_y), dtype = np.float64)
        for j in range(size_y):
            r_y[j] = y_next[j] - y_cur[j] - delta_t * sum_i[j]
        sol[i].set_f_b(r_y)
    r_bc = np.zeros((size_y + size_p), dtype = np.float64)
    y0_tmp = np.zeros((size_y), dtype = np.float64)
    yN_tmp = np.zeros((size_y), dtype = np.float64)
    for j in range(size_y):
        y0_tmp[j] = y[j][0] 
        yN_tmp[j] = y[j][N - 1] 
    bvp_dae._abvp_r(y0_tmp, yN_tmp, p, r_bc)
    sol[N - 1].set_f_N(r_bc)
    start_index_F = (N - 1) * (size_y + m * (size_y + size_z))
    for j in range(size_y + size_p):
        F[start_index_F + j] = r_bc[j]
    for i in range(N - 1):
        start_index_f_a = i * (size_y + m * (size_y + size_z))
        for j in range((m * (size_y + size_z))):
            F[start_index_f_a + j] = sol[i].f_a[j]
        start_index_f_b = start_index_f_a + m * (size_y + size_z)
        for j in range(size_y):
            F[start_index_f_b + j] = sol[i].f_b[j]
    return F, sol

def f_d_jacobian(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan, q, alpha):
    x = np.copy(q)
    size = x.shape[0]
    e = np.eye(size)
    F_s_x, _ = F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan, x, alpha)
    M = np.zeros((size, size), dtype = np.float64)
    eps = np.spacing(1)
    for i in range(size):
        h = math.sqrt(eps) * (max(1, abs(x[i])))
        x_new = x + h * e[:, i]
        F_s_x_new, _ = F_q(bvp_dae, size_y, size_z, size_p, m, N, rk, tspan, x_new, alpha)
        for j in range(size):
            M[j][i] = (F_s_x_new[j] - F_s_x[j]) / h
    return M

def Jacobian_construct(bvp_dae, size_y, size_z, size_p, m, N, rk, alpha, sol):

    b = rk.b
    a = rk.A
    p0 = sol[N - 1].p

    for i in range(N - 1):
        for j in range(m):
            y = sol[i].get_y_tilda(j)
            z = sol[i].get_z_tilda(j)
            Dh = np.zeros((size_y, (size_y + size_z + size_p)), dtype = np.float64)
            Dg = np.zeros((size_z, (size_y + size_z + size_p)), dtype = np.float64)
            bvp_dae._abvp_Df(y, z, p0, Dh)
            bvp_dae._abvp_Dg(y, z, p0, alpha, Dg)
            sol[i].set_Jacobian(a, b[j], Dh, Dg, j)
    Dr = Dh = np.zeros(((size_y + size_p), (size_y + size_y + size_p)), dtype = np.float64)
    y0 = sol[0].get_y()
    yM = sol[0].get_y()
    bvp_dae._abvp_Dr(y0, yM, p0, Dr)
    r_y0 = np.zeros((size_y + size_p, size_y), dtype = np.float64)
    r_yM = np.zeros((size_y + size_p, size_y), dtype = np.float64)
    r_p0 = np.zeros((size_y + size_p, size_p), dtype = np.float64)
    for i in range(size_y + size_p):
        for j in range(size_y):
            r_y0[i][j] = Dr[i][j]
        for j in range(size_y):
            r_yM[i][j] = Dr[i][j + size_y]
        for j in range(size_y):
            r_p0[i][j] = Dr[i][j + (size_y + size_y)]
    sol[0].set_B(r_y0)
    sol[N - 1].set_B(r_yM)
    sol[N - 1].set_VN(r_p0)

    for i in range(N - 1):
        sol[i].update_Jacobian()
    sol[N - 1].set_HN(sol[N - 1].V_N)
    sol[N - 1].set_bN(sol[N - 1].f_N)

def qr_decomposition(size_s, size_p, N, sol):
    sol[0].C_tilda = sol[0].C
    sol[0].G_tilda = sol[0].A
    sol[0].H_tilda = sol[0].H
    sol[0].b_tilda = sol[0].b
    for i in range(N - 2):
        C_tilda_A = np.concatenate((sol[i].C_tilda, sol[i + 1].A), axis = 0)
        Q, R = np.linalg.qr(C_tilda_A, mode = 'complete')
        sol[i].R = R[0 : size_s, :]
        zero_C = np.concatenate((np.zeros((size_s, size_s), dtype = np.float64), sol[i + 1].C), axis = 0)
        EC = np.dot(Q.T, zero_C)
        sol[i].E = EC[0 : size_s, 0 : size_s]
        sol[i + 1].C_tilda = EC[size_s : 2 * size_s, 0 : size_s]
        G_tilda_zero = np.concatenate((sol[i].G_tilda, np.zeros((size_s, size_s), dtype = np.float64)), axis = 0)
        GG = np.dot(Q.T, G_tilda_zero)
        sol[i].G = GG[0 : size_s, 0 : size_s]
        sol[i + 1].G_tilda = GG[size_s : 2 * size_s, 0 : size_s]
        H_tilda_H = np.concatenate((sol[i].H_tilda, sol[i + 1].H), axis = 0)
        JH = np.dot(Q.T, H_tilda_H)
        sol[i].K = JH[0 : size_s, 0 : size_p]
        sol[i + 1].H_tilda = JH[size_s : 2 * size_s, 0 : size_p]
        b_tilda_b = np.concatenate((sol[i].b_tilda, sol[i + 1].b), axis = 0)
        db = np.dot(Q.T, b_tilda_b)
        sol[i].d = db[0 : size_s]
        sol[i + 1].b_tilda = db[size_s : 2 * size_s]
    final_block_up = np.concatenate((sol[N - 2].C_tilda, sol[N - 2].G_tilda, sol[N - 2].H_tilda), axis = 1)
    H_N = np.asarray(sol[N - 1].H_N)
    final_block_down = np.concatenate((sol[N - 1].B, sol[0].B, H_N), axis = 1)
    final_block = np.concatenate((final_block_up, final_block_down), axis = 0)
    Q, R = np.linalg.qr(final_block, mode = 'complete')
    sol[N - 2].R = R[0 : size_s, 0 : size_s]
    sol[N - 2].G = R[0 : size_s, size_s : 2 * size_s]
    sol[N - 2].K = R[0 : size_s, 2 * size_s : 2 * size_s + size_p]
    sol[N - 1].R = R[size_s : 2 * size_s, size_s : 2 * size_s]
    sol[N - 1].K = R[size_s : 2 * size_s, 2 * size_s : 2 * size_s + size_p]
    sol[N - 1].Rp = R[2 * size_s : 2 * size_s + size_p, 2 * size_s : 2 * size_s + size_p]

    b_N = np.asarray(sol[N - 1].b_N)
    b_tilda_b = np.concatenate((sol[N - 2].b_tilda, b_N), axis = 0)
    d = np.dot(Q.T, b_tilda_b)
    sol[N - 2].d = d[0 : size_s]
    sol[N - 1].d = d[size_s : 2 * size_s]
    sol[N - 1].dp = d[2 * size_s : 2 *size_s + size_p]

def backward_substitution(N, sol):
    delta_p = np.linalg.solve(sol[N - 1].Rp, sol[N - 1].dp)
    sol[N - 1].delta_p = delta_p
    delta_y1 = np.linalg.solve(sol[N - 1].R, (sol[N - 1].d - np.dot(sol[N - 1].K, delta_p)))
    sol[0].set_delta_y(delta_y1)
    b_yN = sol[N - 2].d - np.dot(sol[N - 2].G, delta_y1) - np.dot(sol[N - 2].K, delta_p)
    delta_yN = np.linalg.solve(sol[N - 2].R, b_yN)
    sol[N - 1].set_delta_y(delta_yN)
    for i in range(N - 2, 0, -1):
        b_yi = sol[i - 1].d - np.dot(sol[i - 1].E, sol[i + 1].delta_y) - np.dot(sol[i - 1].G, delta_y1) - np.dot(sol[i - 1].K, delta_p)
        delta_yi = np.linalg.solve(sol[i - 1].R, b_yi)
        sol[i].set_delta_y(delta_yi)

    for i in range(N - 1):
        W_inv = np.linalg.inv(sol[i].W)
        tmp = -sol[i].f_a - np.dot(sol[i].J, sol[i].delta_y) - np.dot(sol[i].V, delta_p)
        delta_k = np.dot(W_inv, tmp)
        sol[i].set_delta_k(delta_k)

def get_delta_q(size_y, size_z, size_p, m, N, sol):
    delta_q = np.zeros((N * size_y + (N - 1) * m * (size_y + size_z) + size_p), dtype = np.float64)
    for i in range(N - 1):
        start_index = i * (size_y + m * (size_y + size_z))
        delta_i = np.concatenate((sol[i].delta_y, sol[i].delta_k), axis = 0)
        for j in range ((size_y + m * (size_y + size_z))):
            delta_q[start_index + j] = delta_i[j]
    delta_N = np.concatenate((sol[N - 1].delta_y, sol[N - 1].delta_p), axis = 0)
    start_index = (N - 1) * (size_y + m * (size_y + size_z))
    for j in range(size_y + size_p):
        delta_q[start_index + j] = delta_N[j]
    return delta_q