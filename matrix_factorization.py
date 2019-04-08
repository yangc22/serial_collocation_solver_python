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
        norm_x = math.sqrt(norm_x)
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
