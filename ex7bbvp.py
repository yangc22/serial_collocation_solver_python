# Created by OCP.py
from math import *
from BVPDAEReadWriteData import BVPDAEReadData, BVPDAEWriteData
import numpy as np

class bvp_dae:

	def __init__(self):
		self.size_y = 4
		self.size_z = 1
		self.size_p = 4
		self.size_inequality = 0
		self.tolerance = 1e-06
		self.maximum_nodes = 1000
		self.maximum_newton_iterations = 200
		self.maximum_mesh_refinements = 10
		error, t0, y0, z0, p0 = BVPDAEReadData("ex7a.data")
		if (error != 0):
			print("Unable to read input file!")
		self.N = t0.shape[0]
		self.T0 = t0
		self.Y0 = None
		if (self.size_y > 0):
			self.Y0 = np.ones((self.N, self.size_y), dtype = np.float64)
		self.Z0 = None
		if (self.size_z > 0):
			self.P0 = np.ones((self.N, self.size_z), dtype = np.float64)
		self.P0 = None
		if (self.size_p > 0):
			self.P0 = np.ones((self.size_p), dtype = np.float64)
		_pack_YZP(self.Y0, self.Z0, self.P0, y0, z0, p0);

	def _pack_YZP(self, _Y, _Z, _P, y0, z0, p0):
		if (_Y != None) and (y0 != None):
			_n = self.N
			_m = y0.shape[1] if y0.shape[1] < _Y.shape[1] else _Y.shape[1]
			for i in range(_n):
				for j in range(_m):
					_Y[i][j] = y0[i][j]
		if (_Z != None) and (z0 != None):
			_n = self.N
			# only read in enough dtat to fill the controls
			_m = z0.shape[1] if z0.shape[1] < 1 else 1
			for i in range(_n):
				for j in range(_m):
					_Z[i][j] = z0[i][j]
		if (_P != None) and (p0 != None):
			_n = p0.shape[0] if p0.shape[0] < _P.shape[0] else _P.shape[0]
			for i in range(_n):
				_P[i] = p0[i]

	def _abvp_f(self, _y, _z, _p, _f):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_f[0] = -u1*(x1 + 0.25) - 2.0*x1 + (x2 + 0.5)*exp(25.0*x1/(x1 + 2.0)) - 0.5
		_f[1] = -x2 - (x2 + 0.5)*exp(25.0*x1/(x1 + 2.0)) + 0.5
		_f[2] = -_lambda1*(-u1 + (x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) - 2.0) + _lambda2*(x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) - 1.0*x1
		_f[3] = -_lambda1*exp(25.0*x1/(x1 + 2.0)) - _lambda2*(-exp(25.0*x1/(x1 + 2.0)) - 1) - 1.0*x2


	def _abvp_g(self, _y, _z, _p, _alpha, _g):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_g[0] = _lambda1*(-x1 - 0.25) + 0.1*u1


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
		_r[0] = x1 - 0.05
		_r[1] = x2
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
		_Df[0][0] = -u1 + (x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) - 2.0
		_Df[0][1] = exp(25.0*x1/(x1 + 2.0))
		_Df[0][4] = -x1 - 0.25
		_Df[1][0] = -(x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0))
		_Df[1][1] = -exp(25.0*x1/(x1 + 2.0)) - 1
		_Df[2][0] = -_lambda1*((x2 + 0.5)*(50.0*x1/pow(x1 + 2.0, 3) - 50.0/pow(x1 + 2.0, 2))*exp(25.0*x1/(x1 + 2.0)) + (x2 + 0.5)*pow(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0), 2)*exp(25.0*x1/(x1 + 2.0))) + _lambda2*(x2 + 0.5)*(50.0*x1/pow(x1 + 2.0, 3) - 50.0/pow(x1 + 2.0, 2))*exp(25.0*x1/(x1 + 2.0)) + _lambda2*(x2 + 0.5)*pow(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0), 2)*exp(25.0*x1/(x1 + 2.0)) - 1.0
		_Df[2][1] = -_lambda1*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) + _lambda2*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0))
		_Df[2][2] = u1 - (x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) + 2.0
		_Df[2][3] = (x2 + 0.5)*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0))
		_Df[2][4] = _lambda1
		_Df[3][0] = -_lambda1*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0)) + _lambda2*(-25.0*x1/pow(x1 + 2.0, 2) + 25.0/(x1 + 2.0))*exp(25.0*x1/(x1 + 2.0))
		_Df[3][1] = -1.00000000000000
		_Df[3][2] = -exp(25.0*x1/(x1 + 2.0))
		_Df[3][3] = exp(25.0*x1/(x1 + 2.0)) + 1


	def _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):
		x1 = _y[0]
		x2 = _y[1]
		_lambda1 = _y[2]
		_lambda2 = _y[3]
		u1 = _z[0]
		_Dg[0][0] = -_lambda1
		_Dg[0][2] = -x1 - 0.25
		_Dg[0][4] = 0.100000000000000


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

# Kirk, D. E., Optimal Control Theory: An Introduction,
# Prentice-Hall, 1970.

StateVariables 			= [x1, x2];
ControlVariables 		= [u1];
InitialConstraints 		= [x1-0.05, x2];
TerminalConstraints		= [x1, x2];
CostFunctional 			= 0.5*(x1*x1 + x2*x2 + 0.1*u1*u1);

DifferentialEquations 	= [-2.0*(x1 + 0.25) + (x2 + 0.5)*exp(25.0*x1/(x1 + 2.0)) - (x1 + 0.25)*u1, 0.5 - x2 - (x2 + 0.5)*exp(25.0*x1/(x1 + 2.0))];

#StateEstimate			= [0.05 - 0.09*t/0.78, -t*(0.78-t)];
#ControlEstimate			= [1 - t/0.78];

InitialTime				= 0.0;
FinalTime				= 0.78;
Nodes					= 79;

Tolerance				= 1.0e-6;
InputFile				= "ex7a.data";
OutputFile				= "ex7b.data";
Display					= NO;


'''

