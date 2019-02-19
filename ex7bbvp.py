# Created by OCP.py
from math import *
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
		MSHOOTDAEData _data = MSHOOTDAEReadData("ex7a.data");
		if (_data.error != 0) { RuntimeWarning("Unable to read input file"); ABVPDAEDelete(_bvp); MSHOOTDAEDelete(_m); return _data.error; }
		Vector _T0 = _data.T;
		Matrix _Y0 = NULL; if (_ny > 0) { _Y0 = MatrixNew(_T0->r, _ny); MatrixSetAllTo(_Y0, 1.0); }
		Matrix _Z0 = NULL; if (_nz > 0) { _Z0 = MatrixNew(_T0->r, _nz); MatrixSetAllTo(_Z0, 10.0); }
		Vector _P0 = NULL; if (_np > 0) { _P0 = VectorNew(_np); VectorSetAllTo(_P0, 10.0); }
		_pack_YZP(_Y0, _Z0, _P0, _data);
		if (_data.Y != NULL) MatrixDelete(_data.Y);
		if (_data.Z != NULL) MatrixDelete(_data.Z);
		if (_data.P != NULL) VectorDelete(_data.P);

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

