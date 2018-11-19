# Created by OCP.py
from math import *
import numpy as np

class bvp_dae:

	def __init__(self) :
		self.size_y = 9
		self.size_z = 12
		self.size_p = 9
		self.size_inequality = 5
		self.tolerance = 0.001
		self.maximum_nodes = 5000
		self.maximum_newton_iterations = 500
		self.maximum_mesh_refinements = 50
	MSHOOTDAEData _data = MSHOOTDAEReadData("ex16.data");
	if (_data.error != 0) { RuntimeWarning("Unable to read input file"); ABVPDAEDelete(_bvp); MSHOOTDAEDelete(_m); return _data.error; }
	Vector _T0 = _data.T;
	Matrix _Y0 = NULL; if (_ny > 0) { _Y0 = MatrixNew(_T0->r, _ny); MatrixSetAllTo(_Y0, 1.0); }
	Matrix _Z0 = NULL; if (_nz > 0) { _Z0 = MatrixNew(_T0->r, _nz); MatrixSetAllTo(_Z0, 10.0); }
	Vector _P0 = NULL; if (_np > 0) { _P0 = VectorNew(_np); VectorSetAllTo(_P0, 10.0); }
	_pack_YZP(_Y0, _Z0, _P0, _data);
	if (_data.Y != NULL) MatrixDelete(_data.Y);
	if (_data.Z != NULL) MatrixDelete(_data.Z);
	if (_data.P != NULL) VectorDelete(_data.P);
	_solution_estimate(_T0, _Y0, _Z0, _P0);

	def _solution_estimate(self, _T, _Y, _Z, _P):
		N = _T.shape[0]
		for i in range(N):
			t = _T[i];

		if (_P.shape[0] != 0):
			for i in range(self.size_p):
				_P[0] = 1.0;

	def _abvp_f(self, _y, _z, _p, _f):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		x1 = _y[0]
		x2 = _y[1]
		x3 = _y[2]
		x4 = _y[3]
		_lambda1 = _y[4]
		_lambda2 = _y[5]
		_lambda3 = _y[6]
		_lambda4 = _y[7]
		_gamma1 = _y[8]
		u1 = _z[0]
		u2 = _z[1]
		_mu1 = _z[2]
		_mu2 = _z[3]
		_mu3 = _z[4]
		_mu4 = _z[5]
		_mu5 = _z[6]
		_nu1 = _z[7]
		_nu2 = _z[8]
		_nu3 = _z[9]
		_nu4 = _z[10]
		_nu5 = _z[11]
		p = _p[0]
		_f[0] = p*x3
		_f[1] = p*x4
		_f[2] = 1.0*p*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_f[3] = 1.0*p*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_f[4] = -2*p*x1
		_f[5] = -1.0*_lambda3*p*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 2*p*x2
		_f[6] = -_lambda1*p - 1.0*_lambda3*p*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*p*x3
		_f[7] = -_lambda2*p - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*p*x4
		_f[8] = _lambda1*x3 + _lambda2*x4 + 1.0*_lambda3*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - _mu5 + rho + pow(u1, 2) + pow(u2, 2) + pow(x1, 2) + pow(x2, 2) + pow(x3, 2) + pow(x4, 2)


	def _abvp_g(self, _y, _z, _p, _alpha, _g):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		x1 = _y[0]
		x2 = _y[1]
		x3 = _y[2]
		x4 = _y[3]
		_lambda1 = _y[4]
		_lambda2 = _y[5]
		_lambda3 = _y[6]
		_lambda4 = _y[7]
		_gamma1 = _y[8]
		u1 = _z[0]
		u2 = _z[1]
		_mu1 = _z[2]
		_mu2 = _z[3]
		_mu3 = _z[4]
		_mu4 = _z[5]
		_mu5 = _z[6]
		_nu1 = _z[7]
		_nu2 = _z[8]
		_nu3 = _z[9]
		_nu4 = _z[10]
		_nu5 = _z[11]
		p = _p[0]
		_g[0] = 1.0*_lambda3*p*(Eye_2 + pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + _mu1 - _mu2 + 2*p*u1
		_g[1] = 1.0*_lambda3*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*p*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + _mu3 - _mu4 + 2*p*u2
		_g[2] = _nu1 + u1 - 10.0
		_g[3] = _nu2 - u1 - 10.0
		_g[4] = _nu3 + u2 - 10.0
		_g[5] = _nu4 - u2 - 10.0
		_g[6] = _nu5 - p
		_g[7] = _mu1 + _nu1 - sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2))
		_g[8] = _mu2 + _nu2 - sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2))
		_g[9] = _mu3 + _nu3 - sqrt(2*_alpha + pow(_mu3, 2) + pow(_nu3, 2))
		_g[10] = _mu4 + _nu4 - sqrt(2*_alpha + pow(_mu4, 2) + pow(_nu4, 2))
		_g[11] = _mu5 + _nu5 - sqrt(2*_alpha + pow(_mu5, 2) + pow(_nu5, 2))


	def _abvp_r(self, _y0, _y1, _p, _r):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		p = _p[0]
		_kappa_i1 = _p[1]
		_kappa_i2 = _p[2]
		_kappa_i3 = _p[3]
		_kappa_i4 = _p[4]
		_kappa_f1 = _p[5]
		_kappa_f2 = _p[6]
		_kappa_f3 = _p[7]
		_kappa_f4 = _p[8]
		x1 = _y0[0]
		x2 = _y0[1]
		x3 = _y0[2]
		x4 = _y0[3]
		_lambda1 = _y0[4]
		_lambda2 = _y0[5]
		_lambda3 = _y0[6]
		_lambda4 = _y0[7]
		_gamma1 = _y0[8]
		# initial conditions
		_r[0] = -THETA_10 + x1
		_r[1] = -THETA_20 + x2
		_r[2] = x3
		_r[3] = x4
		_r[4] = _kappa_i1 + _lambda1
		_r[5] = _kappa_i2 + _lambda2
		_r[6] = _kappa_i3 + _lambda3
		_r[7] = _kappa_i4 + _lambda4
		_r[8] = _gamma1
		# final conditions
		x1 = _y1[0]
		x2 = _y1[1]
		x3 = _y1[2]
		x4 = _y1[3]
		_lambda1 = _y1[4]
		_lambda2 = _y1[5]
		_lambda3 = _y1[6]
		_lambda4 = _y1[7]
		_gamma1 = _y1[8]
		_r[9] = -THETA_1f + x1
		_r[10] = -THETA_2f + x2
		_r[11] = x3
		_r[12] = x4
		_r[13] = -_kappa_f1 + _lambda1
		_r[14] = -_kappa_f2 + _lambda2
		_r[15] = -_kappa_f3 + _lambda3
		_r[16] = -_kappa_f4 + _lambda4
		_r[17] = _gamma1


	def _abvp_Df(self, _y, _z, _p, _Df):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		x1 = _y[0]
		x2 = _y[1]
		x3 = _y[2]
		x4 = _y[3]
		_lambda1 = _y[4]
		_lambda2 = _y[5]
		_lambda3 = _y[6]
		_lambda4 = _y[7]
		_gamma1 = _y[8]
		u1 = _z[0]
		u2 = _z[1]
		_mu1 = _z[2]
		_mu2 = _z[3]
		_mu3 = _z[4]
		_mu4 = _z[5]
		_mu5 = _z[6]
		_nu1 = _z[7]
		_nu2 = _z[8]
		_nu3 = _z[9]
		_nu4 = _z[10]
		_nu5 = _z[11]
		p = _p[0]
		_Df[0][2] = p
		_Df[0][21] = x3
		_Df[1][3] = p
		_Df[1][21] = x4
		_Df[2][1] = 1.0*p*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[2][2] = 1.0*p*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[2][3] = 1.0*p*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[2][9] = 1.0*p*(Eye_2 + pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[2][10] = 1.0*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[2][21] = 1.0*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[3][1] = 1.0*p*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[3][2] = 1.0*p*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[3][3] = 1.0*p*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[3][9] = 1.0*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[3][10] = 1.0*p*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[3][21] = 1.0*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[4][0] = -2*p
		_Df[4][21] = -2*x1
		_Df[5][1] = -1.0*_lambda3*p*(-3*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*sin(x2)*cos(x2) + L_1*el_2*m_2*u2*cos(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda3*p*(2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(sin(x2), 2) + 2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(4.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 4*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 3) - 1.0*_lambda4*p*(6.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*sin(x2)*cos(x2) + 6.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*sin(x2)*cos(x2) + 3*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*sin(x2)*cos(x2) + L_1*el_2*m_2*u1*cos(x2) - 2.0*L_1*el_2*m_2*u2*cos(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(sin(x2), 2) + 2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(4.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 4*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 3) - 2*p
		_Df[5][2] = -1.0*_lambda3*p*(-2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(4.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x4*pow(sin(x2), 2) - 2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][3] = -1.0*_lambda3*p*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x4*pow(sin(x2), 2) - 2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][6] = -1.0*p*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][7] = -1.0*p*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][9] = -1.0*L_1*_lambda4*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(Eye_2 + pow(el_2, 2)*m_2)*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][10] = -1.0*L_1*_lambda3*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2.0*L_1*_lambda4*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[5][21] = -1.0*_lambda3*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 2*x2
		_Df[6][1] = -1.0*_lambda3*p*(-2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(4.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x4*pow(sin(x2), 2) - 2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[6][2] = -2.0*L_1*_lambda3*el_2*m_2*p*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2.0*L_1*_lambda4*el_2*m_2*p*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*p
		_Df[6][3] = -2.0*L_1*_lambda3*el_2*m_2*p*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2.0*L_1*_lambda4*el_2*m_2*p*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[6][4] = -p
		_Df[6][6] = -1.0*p*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[6][7] = -1.0*p*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[6][21] = -_lambda1 - 1.0*_lambda3*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*x3
		_Df[7][1] = -1.0*_lambda3*p*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) - 1.0*_lambda4*p*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*pow(sin(x2), 2) + 2*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x4*pow(sin(x2), 2) - 2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Df[7][2] = -2.0*L_1*_lambda3*el_2*m_2*p*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2.0*L_1*_lambda4*el_2*m_2*p*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[7][3] = -2.0*L_1*_lambda3*el_2*m_2*p*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2.0*L_1*_lambda4*el_2*m_2*p*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*p
		_Df[7][5] = -p
		_Df[7][6] = -1.0*p*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[7][7] = -1.0*p*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[7][21] = -_lambda2 - 1.0*_lambda3*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 1.0*_lambda4*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2*x4
		_Df[8][0] = 2*x1
		_Df[8][1] = 1.0*_lambda3*(-pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u2*sin(x2) + L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda3*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) + 1.0*_lambda4*(2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x3, 2)*pow(sin(x2), 2) + 2.0*pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*x3*x4*pow(sin(x2), 2) + pow(L_1, 2)*pow(el_2, 2)*pow(m_2, 2)*pow(x4, 2)*pow(sin(x2), 2) + L_1*el_2*m_2*u1*sin(x2) - 2.0*L_1*el_2*m_2*u2*sin(x2) - L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*cos(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*cos(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) + 2*x2
		_Df[8][2] = _lambda1 + 1.0*_lambda3*(2*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(-2*L_1*el_2*m_2*x3*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*x3
		_Df[8][3] = _lambda2 + 1.0*_lambda3*(2.0*L_1*el_2*m_2*x3*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + 2*L_1*el_2*m_2*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(-2.0*L_1*el_2*m_2*x3*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*x4
		_Df[8][4] = x3
		_Df[8][5] = x4
		_Df[8][6] = 1.0*(L_1*el_2*m_2*pow(x3, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) + 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) + u1*(Eye_2 + pow(el_2, 2)*m_2) - u2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[8][7] = 1.0*(-L_1*el_2*m_2*pow(x3, 2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))*sin(x2) - 2.0*L_1*el_2*m_2*x3*x4*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - L_1*el_2*m_2*pow(x4, 2)*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2) - u1*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2) + u2*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Df[8][9] = 1.0*_lambda3*(Eye_2 + pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*u1
		_Df[8][10] = 1.0*_lambda3*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*u2
		_Df[8][15] = -1


	def _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		x1 = _y[0]
		x2 = _y[1]
		x3 = _y[2]
		x4 = _y[3]
		_lambda1 = _y[4]
		_lambda2 = _y[5]
		_lambda3 = _y[6]
		_lambda4 = _y[7]
		_gamma1 = _y[8]
		u1 = _z[0]
		u2 = _z[1]
		_mu1 = _z[2]
		_mu2 = _z[3]
		_mu3 = _z[4]
		_mu4 = _z[5]
		_mu5 = _z[6]
		_nu1 = _z[7]
		_nu2 = _z[8]
		_nu3 = _z[9]
		_nu4 = _z[10]
		_nu5 = _z[11]
		p = _p[0]
		_Dg[0][1] = 1.0*L_1*_lambda4*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda3*p*(Eye_2 + pow(el_2, 2)*m_2)*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) + 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Dg[0][6] = 1.0*p*(Eye_2 + pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Dg[0][7] = 1.0*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Dg[0][9] = 2*p
		_Dg[0][11] = 1
		_Dg[0][12] = -1
		_Dg[0][21] = 1.0*_lambda3*(Eye_2 + pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*u1
		_Dg[1][1] = 1.0*L_1*_lambda3*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) - 2.0*L_1*_lambda4*el_2*m_2*p*sin(x2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda3*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2) + 1.0*_lambda4*p*(2.0*L_1*el_2*m_2*(Eye_2 + pow(el_2, 2)*m_2)*sin(x2) - 2*L_1*el_2*m_2*(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2)*sin(x2))*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/pow((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2), 2)
		_Dg[1][6] = 1.0*p*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Dg[1][7] = 1.0*p*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2))
		_Dg[1][10] = 2*p
		_Dg[1][13] = 1
		_Dg[1][14] = -1
		_Dg[1][21] = 1.0*_lambda3*(-Eye_2 - L_1*el_2*m_2*cos(x2) - pow(el_2, 2)*m_2)/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 1.0*_lambda4*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2)))/((Eye_2 + pow(el_2, 2)*m_2)*(Eye_1 + Eye_2 + pow(el_1, 2)*m_1 + m_2*(pow(L_1, 2) + 2.0*L_1*el_2*cos(x2) + pow(el_2, 2))) - pow(Eye_2 + L_1*el_2*m_2*cos(x2) + pow(el_2, 2)*m_2, 2)) + 2*u2
		_Dg[2][9] = 1
		_Dg[2][16] = 1
		_Dg[3][9] = -1
		_Dg[3][17] = 1
		_Dg[4][10] = 1
		_Dg[4][18] = 1
		_Dg[5][10] = -1
		_Dg[5][19] = 1
		_Dg[6][20] = 1
		_Dg[6][21] = -1
		_Dg[7][11] = -_mu1/sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2)) + 1
		_Dg[7][16] = -_nu1/sqrt(2*_alpha + pow(_mu1, 2) + pow(_nu1, 2)) + 1
		_Dg[8][12] = -_mu2/sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2)) + 1
		_Dg[8][17] = -_nu2/sqrt(2*_alpha + pow(_mu2, 2) + pow(_nu2, 2)) + 1
		_Dg[9][13] = -_mu3/sqrt(2*_alpha + pow(_mu3, 2) + pow(_nu3, 2)) + 1
		_Dg[9][18] = -_nu3/sqrt(2*_alpha + pow(_mu3, 2) + pow(_nu3, 2)) + 1
		_Dg[10][14] = -_mu4/sqrt(2*_alpha + pow(_mu4, 2) + pow(_nu4, 2)) + 1
		_Dg[10][19] = -_nu4/sqrt(2*_alpha + pow(_mu4, 2) + pow(_nu4, 2)) + 1
		_Dg[11][15] = -_mu5/sqrt(2*_alpha + pow(_mu5, 2) + pow(_nu5, 2)) + 1
		_Dg[11][20] = -_nu5/sqrt(2*_alpha + pow(_mu5, 2) + pow(_nu5, 2)) + 1


	def _abvp_Dr(self, _y0, _y1, _p, _Dr):
		rho =1.0e3
		THETA_10 =0.0
		THETA_20 =-2.0
		THETA_1f =1.0
		THETA_2f =-1.0
		L_1 =0.4
		L_2 =0.4
		m_1 =0.5
		m_2 =0.5
		Eye_1 =0.1
		Eye_2 =0.1
		el_1 =0.2
		el_2 =0.2
		p = _p[0]
		_kappa_i1 = _p[1]
		_kappa_i2 = _p[2]
		_kappa_i3 = _p[3]
		_kappa_i4 = _p[4]
		_kappa_f1 = _p[5]
		_kappa_f2 = _p[6]
		_kappa_f3 = _p[7]
		_kappa_f4 = _p[8]
		x1 = _y0[0]
		x2 = _y0[1]
		x3 = _y0[2]
		x4 = _y0[3]
		_lambda1 = _y0[4]
		_lambda2 = _y0[5]
		_lambda3 = _y0[6]
		_lambda4 = _y0[7]
		_gamma1 = _y0[8]
		# initial conditions
		_Dr[0][0] = 1
		_Dr[1][1] = 1
		_Dr[2][2] = 1
		_Dr[3][3] = 1
		_Dr[4][4] = 1
		_Dr[4][19] = 1
		_Dr[5][5] = 1
		_Dr[5][20] = 1
		_Dr[6][6] = 1
		_Dr[6][21] = 1
		_Dr[7][7] = 1
		_Dr[7][22] = 1
		_Dr[8][8] = 1
		# final conditions
		x1 = _y1[0]
		x2 = _y1[1]
		x3 = _y1[2]
		x4 = _y1[3]
		_lambda1 = _y1[4]
		_lambda2 = _y1[5]
		_lambda3 = _y1[6]
		_lambda4 = _y1[7]
		_gamma1 = _y1[8]
		_Dr[9][9] = 1
		_Dr[10][10] = 1
		_Dr[11][11] = 1
		_Dr[12][12] = 1
		_Dr[13][13] = 1
		_Dr[13][23] = -1
		_Dr[14][14] = 1
		_Dr[14][24] = -1
		_Dr[15][15] = 1
		_Dr[15][25] = -1
		_Dr[16][16] = 1
		_Dr[16][26] = -1
		_Dr[17][17] = 1


'''

# Dissanayake, M., Goh, C. J., and Phan-Thien, N.,
# Time-optimal trajectories for robot manipulators, 
# Robotica, Vol. 9, pp. 131--138, 1991.

Constants = [rho = 1.0e3, THETA_10=0.0, THETA_20=-2.0, THETA_1f=1.0, THETA_2f=-1.0,
	L_1=0.4, L_2=0.4, m_1=0.5, m_2=0.5, Eye_1=0.1, Eye_2=0.1, el_1=0.2, el_2=0.2];
StateVariables 			= [x1, x2, x3, x4];
ControlVariables 		= [u1, u2];
ParameterVariables		= [p];
InitialConstraints 		= [x1-THETA_10,	x2-THETA_20, x3, x4];
TerminalConstraints 	= [x1-THETA_1f, x2-THETA_2f, x3, x4];
#TerminalPenalty	= (x1-THETA_1f)*(x1-THETA_1f) + (x2-THETA_2f)*(x2-THETA_2f) + x3*x3 + x4*x4;
CostFunctional 			= p*(rho + x1*x1 + x2*x2 + x3*x3 + x4*x4 + u1*u1 + u2*u2);
cs1 = cos(x2);
H_11 = Eye_1+Eye_2+m_1*el_1*el_1+m_2*(L_1*L_1+el_2*el_2+2.0*L_1*el_2*cs1);
H_12 = Eye_2+m_2*el_2*el_2+m_2*L_1*el_2*cs1;
H_22 = Eye_2+m_2*el_2*el_2;
h = m_2*L_1*el_2*sin(x2);
delta = 1.0/(H_11*H_22-H_12*H_12);
	
DifferentialEquations 	= [p*x3, p*x4, 
	p*delta*(2.0*h*H_22*x3*x4+h*H_22*x4*x4+h*H_12*x3*x3+H_22*u1-H_12*u2),
	p*delta*(-2.0*h*H_12*x3*x4-h*H_11*x3*x3-h*H_12*x4*x4+H_11*u2-H_12*u1)];

ParameterEstimate		= [1.0];

InequalityConstraints	= [u1 - 10.0, -10.0 - u1, u2 - 10.0, -10.0 - u2, -p];

Tolerance				= 1.0e-3;
InputFile				= "ex16.data";
OutputFile				= "ex17.data";
Display					= NO;
MaximumNewtonIterations = 500;
MaximumMeshRefinements	= 50;
MaximumNodes            = 5000;

'''

