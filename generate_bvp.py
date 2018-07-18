from sympy import *

x1, x2 = symbols('x1 x2')
u1 = symbols('u1')
cf = u1 * u1 * 0.5
de = [x2, u1]
ic = [x1 - 1, x2 - 1]
tc = [x1, x2]

size_x = 2
_lambda = symbols('_lambda1 _lambda2')
hamiltonian = cf
for i in range(size_x):
    hamiltonian += _lambda[i] * de[i]

y = [x1 x2 _lambda1 _lambda2]
h = diff(hamiltonian, y[0])