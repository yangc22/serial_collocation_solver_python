def ODE_h(y, z, p, alpha):
    lambda1 = y[3]
    x2 = y[2]
    u1 = z[0]
    h = [x2, u1, 0, -lambda1]
    return h