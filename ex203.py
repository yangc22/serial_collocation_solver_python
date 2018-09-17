from math import *

class bvp_dae:

	def _abvp_f(self, _y, _z, _p, _f):
		gamma =0.5
		x = _y[0]
		y = _y[1]
		theta = _y[2]
		_lambda1 = _y[3]
		_lambda2 = _y[4]
		_lambda3 = _y[5]
		_gamma1 = _y[6]
		u = _z[0]
		p = _p[0]
		_f[0] = (p + 18)*cos(theta)
		_f[1] = (p + 18)*sin(theta)
		_f[2] = u*(p + 18)
		_f[3] = 0
		_f[4] = 0
		_f[5] = _lambda1*(p + 18)*sin(theta) - _lambda2*(p + 18)*cos(theta)
		_f[6] = _lambda1*cos(theta) + _lambda2*sin(theta) + _lambda3*u + gamma*pow(u, 2) + 1


	def _abvp_g(self, _y, _z, _p, _alpha, _g):
		gamma =0.5
		x = _y[0]
		y = _y[1]
		theta = _y[2]
		_lambda1 = _y[3]
		_lambda2 = _y[4]
		_lambda3 = _y[5]
		_gamma1 = _y[6]
		u = _z[0]
		p = _p[0]
		_g[0] = _lambda3*(p + 18) + 2*gamma*u*(p + 18)


	def _abvp_r(self, _y0, _y1, _p, _r):
		gamma =0.5
		p = _p[0]
		_kappa_i1 = _p[1]
		_kappa_i2 = _p[2]
		_kappa_i3 = _p[3]
		_kappa_f1 = _p[4]
		_kappa_f2 = _p[5]
		_kappa_f3 = _p[6]
		x = _y0[0]
		y = _y0[1]
		theta = _y0[2]
		_lambda1 = _y0[3]
		_lambda2 = _y0[4]
		_lambda3 = _y0[5]
		_gamma1 = _y0[6]
		# initial conditions
		_r[0] = x + 9
		_r[1] = y
		_r[2] = theta
		_r[3] = _kappa_i1 + _lambda1
		_r[4] = _kappa_i2 + _lambda2
		_r[5] = _kappa_i3 + _lambda3
		_r[6] = _gamma1
		# final conditions
		x = _y1[0]
		y = _y1[1]
		theta = _y1[2]
		_lambda1 = _y1[3]
		_lambda2 = _y1[4]
		_lambda3 = _y1[5]
		_gamma1 = _y1[6]
		_r[7] = x - 9
		_r[8] = y
		_r[9] = theta
		_r[10] = -_kappa_f1 + _lambda1
		_r[11] = -_kappa_f2 + _lambda2
		_r[12] = -_kappa_f3 + _lambda3
		_r[13] = _gamma1


	def _abvp_Df(self, _y, _z, _p, _Df):
		gamma =0.5
		x = _y[0]
		y = _y[1]
		theta = _y[2]
		_lambda1 = _y[3]
		_lambda2 = _y[4]
		_lambda3 = _y[5]
		_gamma1 = _y[6]
		u = _z[0]
		p = _p[0]
		_Df[0][2] = -(p + 18)*sin(theta)
		_Df[0][8] = cos(theta)
		_Df[1][2] = (p + 18)*cos(theta)
		_Df[1][8] = sin(theta)
		_Df[2][7] = p + 18
		_Df[2][8] = u
		_Df[5][2] = _lambda1*(p + 18)*cos(theta) + _lambda2*(p + 18)*sin(theta)
		_Df[5][3] = (p + 18)*sin(theta)
		_Df[5][4] = -(p + 18)*cos(theta)
		_Df[5][8] = _lambda1*sin(theta) - _lambda2*cos(theta)
		_Df[6][2] = -_lambda1*sin(theta) + _lambda2*cos(theta)
		_Df[6][3] = cos(theta)
		_Df[6][4] = sin(theta)
		_Df[6][5] = u
		_Df[6][7] = _lambda3 + 2*gamma*u


	def _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):
		gamma =0.5
		x = _y[0]
		y = _y[1]
		theta = _y[2]
		_lambda1 = _y[3]
		_lambda2 = _y[4]
		_lambda3 = _y[5]
		_gamma1 = _y[6]
		u = _z[0]
		p = _p[0]
		_Dg[0][5] = p + 18
		_Dg[0][7] = 2*gamma*(p + 18)
		_Dg[0][8] = _lambda3 + 2*gamma*u


	def _abvp_Dr(self, _y0, _y1, _p, _Dr):
		gamma =0.5
		p = _p[0]
		_kappa_i1 = _p[1]
		_kappa_i2 = _p[2]
		_kappa_i3 = _p[3]
		_kappa_f1 = _p[4]
		_kappa_f2 = _p[5]
		_kappa_f3 = _p[6]
		x = _y0[0]
		y = _y0[1]
		theta = _y0[2]
		_lambda1 = _y0[3]
		_lambda2 = _y0[4]
		_lambda3 = _y0[5]
		_gamma1 = _y0[6]
		# initial conditions
		_Dr[0][0] = 1
		_Dr[1][1] = 1
		_Dr[2][2] = 1
		_Dr[3][3] = 1
		_Dr[3][15] = 1
		_Dr[4][4] = 1
		_Dr[4][16] = 1
		_Dr[5][5] = 1
		_Dr[5][17] = 1
		_Dr[6][6] = 1
		# final conditions
		x = _y1[0]
		y = _y1[1]
		theta = _y1[2]
		_lambda1 = _y1[3]
		_lambda2 = _y1[4]
		_lambda3 = _y1[5]
		_gamma1 = _y1[6]
		_Dr[7][7] = 1
		_Dr[8][8] = 1
		_Dr[9][9] = 1
		_Dr[10][10] = 1
		_Dr[10][18] = -1
		_Dr[11][11] = 1
		_Dr[11][19] = -1
		_Dr[12][12] = 1
		_Dr[12][20] = -1
		_Dr[13][13] = 1


'''

# Homotopy Algorithm for Optimal Control Problems
# with a Second-order State Constraint
# Audrey Hermant
# Appl Math Optim (2010) 61: 85–127

Constants               = [gamma = 0.5];
ParameterVariables      = [p];
StateVariables          = [x, y, theta];
ControlVariables        = [u];

InitialConstraints      = [x + 9, y, theta];
TerminalConstraints     = [x - 9, y, theta];

CostFunctional          = (p + 18)*(1 + gamma*u*u);

DifferentialEquations   = [(p + 18)*cos(theta), (p + 18)*sin(theta), (p + 18)*u];

StateVariableInequalityConstraints = [10*cos(0.2*x) - 0.002*x*x*x*x*sin(x-1.25) - y];

Nodes                   = 101;
Tolerance               = 1.0e-6;

StateEstimate           = [18*t-9, 0, -1.5];
ControlEstimate         = [0.001];
ParameterEstimate       = [0];

OutputFile              = "ex203.data";

'''

