from sympy import *

StateVariables          = "[x1, x2]"
ControlVariables        = "[u1]"
InitialConstraints      = "[x1 - 1.0, x2 - 1.0]"
TerminalConstraints     = "[x1, x2]"
CostFunctional          = "u1*u1*0.5"
DifferentialEquations   = "[x2, u1]"

sv = symbols(StateVariables)
cv = symbols(ControlVariables)
ic = symbols(InitialConstraints)
tc = symbols(TerminalConstraints)
cf = symbols(CostFunctional)
de = symbols(DifferentialEquations)
