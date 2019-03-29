from math import *
import numpy as np

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
    v = np.zeros((n,n), dtype=np.float64)

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
            v[k, i] = xi
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
                s += R[j, i] * v[k, j]
            w[i] = s
        # R(k: n, k: m) += beta * v(k: n) *w(k: m) ^ T
        # a[k: n, k: m] = a[k: n, k: m] + beta * (v * wT)
        for i in range(k, n):
            for j in range(k, m):
                R[i, j] += beta[k] * v[k, i] * w[j]

        for k in range(m - 1, -1, -1):
            # double * v = V->e[k];
            # V is a matrix, v is the element of the kth row
            # double * u = W->e;
            # W is a vector, u is the element of the vector
            # uT = v.transpose() * Q[k: n, k: n]
            for i in range(k, n):
                s = 0.0
                for j in range(k, n):
                    s += Q[j, i] * v[k, j]
                w[i] = s
            # Q[k: n, k: n] = Q[k: n, k: n] + Beta[k] * (v * uT)
            for i in range(k, n):
                for j in range(k, n):
                    Q[i, j] += beta[k] * v[k, i] * w[j]
