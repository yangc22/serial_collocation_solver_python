from sympy import *
'''
def h_ODE(y):
    x1, x2 = symbols('x1 x2')
    u1 = symbols('u1')
    cf = u1 * u1 * 0.5
    de = [x2, u1]
    ic = [x1 - 1, x2 - 1]
    tc = [x1, x2]

    _lambda = symbols('_lambda1 _lambda2')
    hamiltonian = cf

    for i in range(size_x):
        hamiltonian += _lambda[i] * de[i]
    y = []
    y.append(x1)
    y.append(x2)
    y.append(_lambda[0])
    y.append(_lambda[1])
    h = []
    for i in range(len(y)):
        h.append(diff(hamiltonian, y[i]))
    return h

f = lambdify(y, h_ODE(y))
'''

def lagrange_eqs(a):
    x,y,z= symbols('x y z')
    FUNC=x**2-2*x*y**2+z+a*exp(z)
    d_lgrng_1=diff(FUNC,x)
    d_lgrng_2=diff(FUNC,y)
    d_lgrng_3=diff(FUNC,z)
    return [d_lgrng_1,d_lgrng_2,d_lgrng_3]

f = lambdify(((x, y, z),), lagrange_eqs(a))