from math import *

class bvp_dae:

	def _abvp_f(self, _y, _z, _p, _f):
		rho =1.0e-3
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u = _z[0]
		_mu1 = _z[1]
		_mu2 = _z[2]
		_nu1 = _z[3]
		_nu2 = _z[4]
		_f[0] = x2
		_f[1] = u
		_f[2] = 0
		_f[3] = -_lambda1 + 1


	def _abvp_g(self, _y, _z, _p, _alpha, _g):
		rho =1.0e-3
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u = _z[0]
		_mu1 = _z[1]
		_mu2 = _z[2]
		_nu1 = _z[3]
		_nu2 = _z[4]
		_g[0] = _lambda2 + _mu1 - _mu2
		_g[1] = _nu1 + u - 1.0
		_g[2] = _nu2 - u - 1.0
		_g[3] = _mu1 + _nu1 - sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2))
		_g[4] = _mu2 + _nu2 - sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2))


	def _abvp_r(self, _y0, _y1, _p, _r):
		rho =1.0e-3
		_kappa_i1 = _p[0]
		_kappa_i2 = _p[1]
		_kappa_f1 = _p[2]
		x1 = _y0[0]
		x2 = _y0[1]
		_lambda1 = _y0[2]
		_lambda2 = _y0[3]
		# initial conditions
		_r[0] = x1
		_r[1] = x2
		_r[2] = _kappa_i1 + _lambda1
		_r[3] = _kappa_i2 + _lambda2
		# final conditions
		x1 = _y1[0]
		x2 = _y1[1]
		_lambda1 = _y1[2]
		_lambda2 = _y1[3]
		_r[4] = x2
		_r[5] = _lambda1
		_r[6] = -_kappa_f1 + _lambda2


	def _abvp_Df(self, _y, _z, _p, _Df):
		rho =1.0e-3
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u = _z[0]
		_mu1 = _z[1]
		_mu2 = _z[2]
		_nu1 = _z[3]
		_nu2 = _z[4]
		_Df[0][1] = 1
		_Df[1][4] = 1
		_Df[3][2] = -1


	def _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):
		rho =1.0e-3
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u = _z[0]
		_mu1 = _z[1]
		_mu2 = _z[2]
		_nu1 = _z[3]
		_nu2 = _z[4]
		_Dg[0][3] = 1
		_Dg[0][5] = 1
		_Dg[0][6] = -1
		_Dg[1][4] = 1
		_Dg[1][7] = 1
		_Dg[2][4] = -1
		_Dg[2][8] = 1
		_Dg[3][5] = -_mu1/sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2)) + 1
		_Dg[3][7] = -_nu1/sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2)) + 1
		_Dg[4][6] = -_mu2/sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2)) + 1
		_Dg[4][8] = -_nu2/sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2)) + 1


	def _abvp_Dr(self, _y0, _y1, _p, _Dr):
		rho =1.0e-3
		_kappa_i1 = _p[0]
		_kappa_i2 = _p[1]
		_kappa_f1 = _p[2]
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
		_Dr[4][5] = 1
		_Dr[5][6] = 1
		_Dr[6][7] = 1
		_Dr[6][10] = -1


'''

# Zhao, Y., Bryson, A. E. and Slattery, R.,
# Generalized gradient algorithm for trajectory optimization,
# AIAA Journal of Guidance, Control and Dynamics, Vol.  13, pp. 1166--1169, 1990.

Constants               = [rho = 1.0e-3];
StateVariables          = [x1, x2];
ControlVariables        = [u];
InitialConstraints      = [x1, x2];
TerminalConstraints     = [x2];
CostFunctional          = -x2; # + 0.5*rho*u*u;
DifferentialEquations   = [x2, u];
InequalityConstraints   = [u - 1.0, -1.0 - u];
Nodes                   = 101;
Tolerance               = 1.0e-6;
StateEstimate           = [t, t*(1-t)];
ControlEstimate         = [2];
OutputFile              = "ex4.data";
Display                 = NO;
MaximumNodes            = 1000;
MaximumMeshRefinements  = 20;

'''

