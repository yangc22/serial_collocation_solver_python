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

print(hamiltonian)
#y = [x1 x2 _lambda1 _lambda2]
y = []
y.append(x1)
y.append(x2)
y.append(_lambda[0])
y.append(_lambda[1])

h = []
for i in range(len(y)):
    h.append(diff(hamiltonian, y[i]))

print(h)

def lamb(func):
    return lambdify(y, func, modules='numpy')

h = map(lamb, h)

h = h(1)
'''
h_ODE = lambdify(y, h_ODE(hamiltonian, y))

input = [1,2,3,4]
print(h_ODE(input))
'''