import numpy as np
from math import *
import sys


def qr(a):
    """
        QR factorization of self via Householder transformation.
        Returns Q, R such that a = Q*R.
    """
    n, m = a.shape
    if (m > n):
        raise TypeError('qr requires a matrix with columns <= rows')
    Q = np.eye(n)
    V = np.zeros((n, m))
    Beta = np.zeros(m)
    r = m
    # for k in xrange(1, m + 1):
    for k in range(0, m):
        x = np.copy(a[k: n, k])
        if x[0] >= 0:
            s = 1.0
        else:
            s = -1.0
        g = -s * sqrt(x.dot(x))
        v = x
        v[0] = v[0] - g
        vs = v.dot(v)
        if vs < sys.float_info.epsilon:
            r = k
            break
        beta = -2.0 / float(vs)
        wT = v.T.dot(a[k: n, k: m])
        a[k: n, k: m] = a[k: n, k: m] + beta * np.outer(v, wT)  # FIX
        Beta[k] = beta
        V[k: n, k] = np.copy(v)
    # R = CTMatrix(n, m)
    R = np.zeros((n, m))
    R[0: m, 0: m] = np.copy(a[0: m, 0: m])
    # for k in xrange(r, 0, -1):
    for k in range(r - 1, -1, -1):
        v = V[k: n, k]
        uT = v.T.dot(Q[k: n, k: n])
        Q[k: n, k: n] = Q[k: n, k: n] + Beta[k] * np.outer(v, uT)  # FIX
    return Q, R


def qr_cuda(a, cpy):
    '''
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
    '''
    n, m = a.shape
    if m > n:
        raise TypeError('qr requires a matrix with columns <= rows')
    for i in range(n):
        for j in range(m):
            cpy[i, j] = a[i, j]
    q = np.zeros((n, n), dtype=np.float64)
    v = np.zeros((n, m), dtype=np.float64)
    beta = np.zeros(m, dtype=np.float64)
    col = m
    r = np.zeros((n, m), dtype=np.float64)
    vv = np.zeros(n, dtype=np.float64)
    w_t = np.zeros(m, dtype=np.float64)
    u_t = np.zeros(n, dtype=np.float64)
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
            vv[i] = cpy[i, k]
        # x = np.copy(a[k: n, k])
        # if x[0] >= 0:
        #     s = 1.0
        # else:
        #     s = -1.0
        if vv[k] >= 0:
            s = 1.0
        else:
            s = -1.0
        vv_dot = 0
        for i in range(k, n):
            vv_dot += vv[i] * vv[i]
        g = -s * sqrt(vv_dot)
        vv[k] = vv[k] - g
        # vs = vv.dot(vv)
        vs = 0
        for i in range(k, n):
            vs += vv[i] * vv[i]
        if vs < sys.float_info.epsilon:
            col = k
            break
        b = -2.0 / vs
        # w_t = vv[k: n].T.dot(a[k: n, k: m])
        for i in range(k, m):
            w_t[i] = 0
            for j in range(k, n):
                w_t[i] += vv[j] * cpy[j, i]
        # a[k: n, k: m] = a[k: n, k: m] + b * np.outer(vv[k: n], w_t)  # FIX
        for i in range(k, n):
            for j in range(k, m):
                cpy[i, j] += b * vv[i] * w_t[j]
        beta[k] = b
        # v[k: n, k] = np.copy(vv[k: n])
        for i in range(k, n):
            v[i, k] = vv[i]
    # r[0: m, 0: m] = np.copy(a[0: m, 0: m])
    for i in range(m):
        for j in range(m):
            r[i, j] = cpy[i, j]
    # set the other elements as 0s
    for i in range(m, n):
        for j in range(m):
            r[i, j] = 0
    # for k in xrange(r, 0, -1):
    for k in range(col - 1, -1, -1):
        # vv = v[k: n, k]
        for i in range(k, n):
            vv[i] = v[i, k]
        # u_t = vv.T.dot(q[k: n, k: n])
        for i in range(k, n):
            u_t[i] = 0
            for j in range(k, n):
                u_t[i] += vv[j] * q[j, i]
        # q[k: n, k: n] = q[k: n, k: n] + beta[k] * np.outer(vv[k: n], u_t[k: n])  # FIX
        for i in range(k, n):
            for j in range(k, n):
                q[i, j] += beta[k] * vv[i] * u_t[j]
    return q, r

