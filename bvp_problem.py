from math import *

class bvp_dae:

	def _abvp_f(self, _y, _z, _p, _f):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_f[0] = x2
		_f[1] = u1
		_f[2] = 0
		_f[3] = -_lambda1


	def _abvp_g(self, _y, _z, _p, _alpha, _g):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_g[0] = _lambda2 + 1.0*u1


	def _abvp_r(self, _y0, _y1, _p, _r):
		_kappa_i1 = _p[0]
		_kappa_i2 = _p[1]
		_kappa_f1 = _p[2]
		_kappa_f2 = _p[3]
		x1 = _y0[0]
		x2 = _y0[1]
		_lambda1 = _y0[2]
		_lambda2 = _y0[3]
		# initial conditions
		_r[0] = x1 - 1.0
		_r[1] = x2 - 1.0
		_r[2] = _kappa_i1 + _lambda1
		_r[3] = _kappa_i2 + _lambda2
		# final conditions
		x1 = _y1[0]
		x2 = _y1[1]
		_lambda1 = _y1[2]
		_lambda2 = _y1[3]
		_r[4] = x1
		_r[5] = x2
		_r[6] = -_kappa_f1 + _lambda1
		_r[7] = -_kappa_f2 + _lambda2


	def _abvp_Df(self, _y, _z, _p, _Df):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_Df[0][1] = 1
		_Df[1][4] = 1
		_Df[3][2] = -1


	def _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_Dg[0][3] = 1
		_Dg[0][4] = 1.00000000000000


	def _abvp_Dr(self, _y0, _y1, _p, _Dr):
		_kappa_i1 = _p[0]
		_kappa_i2 = _p[1]
		_kappa_f1 = _p[2]
		_kappa_f2 = _p[3]
		x1 = _y0[0]
		x2 = _y0[1]
		_lambda1 = _y0[2]
		_lambda2 = _y0[3]
		# initial conditions
		_Dr[0][0] = 1
		_Dr[1][1] = 1
		_Dr[2][2] = 1
		_Dr[2][8] = 1
		_Dr[3][3] = 1
		_Dr[3][9] = 1
		# final conditions
		x1 = _y1[0]
		x2 = _y1[1]
		_lambda1 = _y1[2]
		_lambda2 = _y1[3]
		_Dr[4][4] = 1
		_Dr[5][5] = 1
		_Dr[6][6] = 1
		_Dr[6][10] = -1
		_Dr[7][7] = 1
		_Dr[7][11] = -1


'''

# Sage, A. P. and White, C. C., Optimum Systems Control,
# Prentice-Hall, 1977.

StateVariables 			= [x1, x2];
ControlVariables 		= [u1];
InitialConstraints 		= [x1 - 1.0, x2 - 1.0];
TerminalConstraints 	= [x1, x2];
CostFunctional 			= u1*u1*0.5;
#TerminalPenalty         = x1*x1;
DifferentialEquations 	= [x2, u1];
Nodes					= 101;
Tolerance				= 1.0e-6;
OutputFile				= "ex2.data";
Display					= YES;

'''

