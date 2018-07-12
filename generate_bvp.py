from sympy import *

StateVariables          = "[x1, x2]"
ControlVariables        = "[u1]"
InitialConstraints      = "[x1 - 1.0, x2 - 1.0]"
TerminalConstraints     = "[x1, x2]"
CostFunctional          = "u1*u1*0.5"
DifferentialEquations   = "[x2, u1]"

# state variables of the OCP
sv = symbols(StateVariables)
# control variables of the OCP
cv = symbols(ControlVariables)
# initial constraints of the OCP
ic = symbols(InitialConstraints)
# terminal constraints of the OCP
tc = symbols(TerminalConstraints)
# cost function of the OCP
cf = symbols(CostFunctional)
# differential equation of the OCP
de = symbols(DifferentialEquations)

# number of the state variables
size_sv = len(sv)
# number of the control variables
size_cv = len(sv)
# number of the initial constraints
size_ic = len(ic)
# number of the terminal constarints
size_tc = len(tc)

# Lagrange multipliers with state variables
_lambda = symarray('_lambda', size_sv)

# Hamiltonian function with Lagrange multipliers
H = cf + _lambda * de

de1 = sv[1]
h = diff(de1, sv[1])
print(de1, sv[1])
print(h)