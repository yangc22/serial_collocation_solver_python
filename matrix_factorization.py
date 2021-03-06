from math import *
import numpy as np
import sys

'''
    QR decomposition of the matrix B, B = Q * R
    Input:
        n : number of rows of B
        m : number of columns of B
        B : matrix to be factorized by QR
            dimension of B : n x m
        Q : holder for the Q matrix, where Q^{T} * Q = I
            dimension of Q : n x n
        R : upper triangular matrix by QR decomposition
            dimension of R : n x m
'''
# n - number of rows
# m - number of columns
# Q - n x n
# R - n x m
# B - n x m
# B = QR


def qr_householder(Q, R, B, n, m, DBL_EPSILON):
    # struct vector * Beta = vector_new(m)
    # struct vector * W = vector_new(n)
    # struct matrix * V = matrix_new(n, n)
    # double * beta = Beta->e
    beta = np.zeros(m, dtype=np.float64)
    w = np.zeros(n, dtype=np.float64)
    v = np.zeros((n, m), dtype=np.float64)

    # Q = I
    for i in range(n):
        for j in range(n):
            if i == j:
                Q[i, j] = 1
            else:
                Q[i, j] = 0
    # R = B
    for i in range(n):
        for k in range(m):
            R[i, k] = B[i, k]
    # find the Householder vectors for each column of R
    for k in range(m):
        # double xi
        norm_x = 0.0
        # double * w = W->e;
        # W is a vector, w is the element of the vector
        # double * v = V->e[k];
        # V is a matrix, v is the element of the kth row
        for i in range(k, n):
            xi = R[i, k]
            norm_x += xi * xi
            v[i, k] = xi
        norm_x = sqrt(norm_x)
        x1 = v[k, k]
        sign_x1 = np.sign(x1)
        # gamma = -1.0 * sign_x1 * norm_x
        g = -1.0 * sign_x1 * norm_x
        v_dot_v = 2.0 * norm_x * (norm_x + sign_x1 * x1)
        if v_dot_v < DBL_EPSILON:
            for i in range(k, m):
                R[i, i] = 0.0
            return
        v[k, k] -= g
        beta[k] = -2.0 / v_dot_v
        # w(k: m) = R(k: n, k: m) ^ T * v(k: n)
        for i in range(k, m):
            s = 0.0
            # FIX: make this row - wise
            for j in range(k, n):
                s += R[j, i] * v[j, k]
            w[i] = s
        # R(k: n, k: m) += beta * v(k: n) *w(k: m) ^ T
        # a[k: n, k: m] = a[k: n, k: m] + beta * (v * wT)
        for i in range(k, n):
            for j in range(k, m):
                R[i, j] += beta[k] * v[i, k] * w[j]
    for k in range(m - 1, -1, -1):
        # double * v = V->e[k];
        # V is a matrix, v is the element of the kth row
        # double * u = W->e;
        # W is a vector, u is the element of the vector
        # uT = v.transpose() * Q[k: n, k: n]
        for i in range(k, n):
            s = 0.0
            for j in range(k, n):
                s += Q[j, i] * v[j, k]
            w[i] = s
        # Q[k: n, k: n] = Q[k: n, k: n] + Beta[k] * (v * uT)
        for i in range(k, n):
            for j in range(k, n):
                Q[i, j] += beta[k] * v[i, k] * w[j]


# def qr(a):
#     m, n = a.shape
#     q = np.eye(m)
#     for i in range(n - (m == n)):
#         h = np.eye(m)
#         v = a[i:, i] / (a[i, i] + np.copysign(np.linalg.norm(a[i:, i]), a[i, i]))
#         v[0] = 1
#         h[i:, i:] = np.eye(a[i:, i].shape[0])
#         h[i:, i:] -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
#         q = np.dot(q, h)
#         a = np.dot(h, a)
#     return q, a

def qr(A):
    m, n = A.shape
    Q = np.eye(m)
    for i in range(n - (m == n)):
        H = np.eye(m)
        H[i:, i:] = make_householder(A[i:, i])
        Q = np.dot(Q, H)
        A = np.dot(H, A)
    return Q, A


def make_householder(a):
    v = a / (a[0] + np.copysign(np.linalg.norm(a), a[0]))
    v[0] = 1
    H = np.eye(a.shape[0])
    H -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
    return H


def qr_cuda(a, cpy, q, r, v, vv, beta, w_t, u_t):
    """
        QR factorization of self via Householder transformation.
        Returns Q, R such that a = Q*R.
        a : matrix, size: n x m
        cpy : matrix, size: n x m
        q : matrix, size: n x n
        v : matrix, size: n x m
        beta : vector, size: m
        r: matrix, size n x m
        vv: vector, size: n(max size)
        w_t: vector, size: m(max size)
        u_t: vector, size: n(max size)
    """
    n, m = a.shape
    if m > n:
        raise TypeError('qr requires a matrix with columns <= rows')
    # copy a matrix
    for i in range(n):
        for j in range(m):
            cpy[i, j] = a[i, j]
    col = m
    # set q as the identity matrix
    for i in range(n):
        for j in range(n):
            if i == j:
                q[i, j] = 1
            else:
                q[i, j] = 0
    # process each column
    for k in range(0, m):
        for i in range(k, n):
            vv[i] = a[i, k]
        if vv[k] >= 0:
            s = 1.0
        else:
            s = -1.0
        vv_dot = 0
        for i in range(k, n):
            vv_dot += vv[i] * vv[i]
        g = -s * sqrt(vv_dot)
        vv[k] = vv[k] - g
        vs = 0
        for i in range(k, n):
            vs += vv[i] * vv[i]
        if vs < sys.float_info.epsilon:
            col = k
            break
        b = -2.0 / vs
        for i in range(k, m):
            w_t[i] = 0
            for j in range(k, n):
                w_t[i] += vv[j] * a[j, i]
        for i in range(k, n):
            for j in range(k, m):
                a[i, j] += b * vv[i] * w_t[j]
        beta[k] = b
        for i in range(k, n):
            v[i, k] = vv[i]
    for i in range(m):
        for j in range(m):
            r[i, j] = a[i, j]
    # set the other elements as 0s
    for i in range(m, n):
        for j in range(m):
            r[i, j] = 0
    for k in range(col - 1, -1, -1):
        for i in range(k, n):
            vv[i] = v[i, k]
        for i in range(k, n):
            u_t[i] = 0
            for j in range(k, n):
                u_t[i] += vv[j] * q[j, i]
        for i in range(k, n):
            for j in range(k, n):
                q[i, j] += beta[k] * vv[i] * u_t[j]
    return q, r


'''
        LU factorization of self.  Returns L, U, P.
        P * a = L * U
'''


def lu(a):
    cpy = np.copy(a)
    n, m = cpy.shape
    if n != m:
        raise TypeError('The input to lu should be a square matrix!')
    P = np.eye(n)
    # L = CTMatrix(n, n)
    # U = CTMatrix(n, n)  # FIX: let A = U,
    L = np.zeros((n, n), dtype=np.float64)
    U = np.zeros((n, n), dtype=np.float64)
    # for j in xrange(n):
    # for each column
    for j in range(n):
        # for i in xrange(j):
        for i in range(j):
            s = 0
            # for k in xrange(i):
            for k in range(i):
                # s = s + L.element[i][k] * U.element[k][j]
                s += L[i, k] * U[k, j]
            U[i, j] = cpy[i, j] - s
        max_u = 0
        row_max = j
        # for i in xrange(j, n):
        for i in range(j, n):
            s = 0
            # for k in xrange(j):
            for k in range(j):
                # s = s + L.element[i][k] * U.element[k][j]
                s += L[i, k] * U[k, j]
            U[i, j] = cpy[i, j] - s
            if i == j:
                max_u = fabs(U[i, i])
            else:
                max_ui = fabs(U[i, j])
                if max_ui > max_u:
                    max_u = max_ui
                    row_max = i
        if row_max != j:  # FIX: this moves too much irrelevant data
            # cpy.swap_rows(j + 1, row_max + 1)  # set U = A
            swap_rows_mat(cpy, j, row_max)
            # L.swap_rows(j + 1, row_max + 1)  # only need to swap columns 1 to j-1
            swap_rows_mat(L, j, row_max)
            # U.swap_rows(j + 1, row_max + 1)  # don't need this if U = A
            swap_rows_mat(U, j, row_max)
            # P.swap_rows(j + 1, row_max + 1)  # store P as a 1d array
            swap_rows_mat(P, j, row_max)
        max_u = U[j, j]
        if fabs(max_u) < 10.0 * sys.float_info.epsilon:
            print('lu(): singular matrix')
        # for i in xrange(j + 1, n):
        for i in range(j + 1, n):
            if fabs(max_u) >= 10.0 * sys.float_info.epsilon:
                L[i, j] = float(U[i, j]) / float(max_u)
            U[i, j] = 0
        L[j, j] = 1
    return L, U, P


'''
    Swap the elements of the ith and jth row of the matrix a.
'''


def swap_rows_mat(a, i, j):
    n, m = a.shape
    if i == j:
        return
    for k in range(m):
        tmp = a[i, k]
        a[i, k] = a[j, k]
        a[j, k] = tmp
    return


'''
    Swap the elements of the ith and jth row of the vector a.
'''


def swap_rows_vec(a, i, j):
    n = a.shape
    if len(n) > 1:
        raise TypeError('The input should be a vector!')
    if i == j:
        return
    tmp = a[i]
    a[i] = a[j]
    a[j] = tmp
    return


'''
    Forward substitution.
    Solve the linear system: Lx = b.
    L is assumed to be lower triangular. (Sanity check)
    b is assumed to be a vector.
'''


def forward_solve_vec(L, b):
    n, m = L.shape
    n_b = b.shape
    min_diag = 1
    for i in range(n):
        min_diag = min(min_diag, abs(L[i, i]))
    if len(n_b) > 1:
        raise TypeError('The input should be a vector in forward_solve_vec!')
    if (n != m) or (n_b[0] != n) or (min_diag < sys.float_info.epsilon):
        raise TypeError('Invalid input to forward_solve!')
    y = np.zeros(n_b[0], dtype=np.float64)
    for i in range(n_b[0]):
        for k in range(i):
            y[i] += L[i, k] * y[k]
        y[i] = float(b[i] - y[i]) / float(L[i, i])
    return y


'''
    Backward substitution.
    Solve the linear system: Ux = y.
    U is assumed to be upper triangular. (Sanity check)
    y is assumed to be a vector.
'''


def backward_solve_vec(U, y):
    n, m = U.shape
    n_y = y.shape
    min_diag = 1
    for i in range(n):
        min_diag = min(min_diag, abs(U[i, i]))
    if len(n_y) > 1:
        raise TypeError('The input should be a vector in backward_solve_vec!')
    if (n != m) or (n_y[0] != n) or (min_diag < sys.float_info.epsilon):
        raise TypeError('Invalid input to backward_solve!')
    x = np.zeros(n_y[0], dtype=np.float64)
    # for i in xrange(y.rows - 1, -1, -1):
    for i in range(n_y[0] - 1, -1, -1):
        for k in range(i, n_y[0]):
            x[i] += U[i, k] * x[k]
        x[i] = float(y[i] - x[i]) / float(U[i, i])
    return x


'''
    Forward substitution.
    Solve the linear system: Lx = b.
    L is assumed to be lower triangular.
    b is assumed to be a matrix.
'''


def forward_solve_mat(L, b):
    n_L, m_L = L. shape
    n_b, m_b = b.shape
    min_diag = 1
    for i in range(n_L):
        min_diag = min(min_diag, abs(L[i, i]))
    if (n_L != m_L) or (n_b != n_L) or (min_diag < sys.float_info.epsilon):
        raise TypeError('Invalid input to forward_solve_mat!')
    y = np.zeros((n_b, m_b), dtype=np.float64)
    # solve each column
    for j in range(m_b):
        for i in range(n_b):
            for k in range(i):
                y[i, j] += L[i, k] * y[k, j]
            y[i, j] = float(b[i, j] - y[i, j]) / float(L[i, i])
    return y


'''
    Backward substitution.
    Solve the linear system: Ux = y.
    U is assumed to be upper triangular.
    y is assumed to be a matrix.
'''


def backward_solve_mat(U, y):
    n_U, m_U = U.shape
    n_y, m_y = y.shape
    min_diag = 1
    for i in range(n_U):
        min_diag = min(min_diag, abs(U[i, i]))
    if (n_U != m_U) or (n_y != n_U) or (min_diag < sys.float_info.epsilon):
        raise TypeError('Invalid input to backward_solve_mat!')
    x = np.zeros((n_y, m_y), dtype=np.float64)
    for j in range(m_y):
        for i in range(n_y - 1, -1, -1):
            for k in range(i, n_y):
                x[i, j] += U[i, k] * x[k, j]
            x[i, j] = float(y[i, j] - x[i, j]) / float(U[i, i])
    return x


'''
    solve the linear system A*x = b
    b is supposed to be a vector
'''


def lu_solve_vec(A, b):
    n_A, m_A = A.shape
    n_b = b.shape
    if len(n_b) > 1:
        raise TypeError('The input to lu_solve_vec should be a vector!')
    if n_A != n_b[0]:
        raise TypeError('Invalid input to lu_solver_vec!')
    L, U, P = lu(A)
    c = np.zeros(n_b[0], dtype=np.float64)
    for i in range(n_A):
        c[i] = 0
        for j in range(m_A):
            c[i] += P[i, j] * b[j]
    y = forward_solve_vec(L, c)
    x = backward_solve_vec(U, y)
    return x


'''
    solve the linear system A*x = b
    b is supposed to be a matrix
'''


def lu_solve_mat(A, b):
    n_A, m_A = A.shape
    n_b, m_b = b.shape
    if n_A != n_b:
        raise TypeError('Invalid input to lu_solver_vec!')
    L, U, P = lu(A)
    c = np.zeros((n_A, m_b), dtype=np.float64)
    for i in range(n_A):
        for j in range(m_A):
            c[i, j] = 0
            for k in range(n_b):
                c[i, j] += P[i, k] * b[k, j]
    y = forward_solve_mat(L, c)
    x = backward_solve_mat(U, y)
    return x
