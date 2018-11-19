"""This file implements the OCP class."""
# TODO: remove unused variables in the C code
import sys
from sympy import *
from sympy import ccode

_raw_line_number_ = 0

class OCP(object):
    """The OCP class"""
    def __init__(self):
        self.x = []
        self.u = []
        self.p = []
        self.phi = 0
        self.L = 0
        self.f = []
        self.c = []
        self.d = []
        self.Gamma = []
        self.Psi = []
        self.t_i = 0
        self.t_f = 1
        self.nodes = 100
        self.input_file = ""
        self.output_file = "ocp.data"
        self.tolerance = 1.0e-6
        self.maximum_nodes = 1000
        self.maximum_newton_iterations = 200
        self.method = "mshootdae"
        self.constants = None
        self.constants_value = None
        self.display = 0
        self.maximum_mesh_refinements = 10
        self.state_estimate = None
        self.control_estimate = None
        self.parameter_estimate = None

    def make_abvp(self):
        '''
            Define the variables for OCP problem.
            y : ODE variables, [x lambda gamma]
            z : DAE variables, [u theta mu nu]
                theta : variables for equality constraints
                mu : lagrange multipliers for inequality constraints
                nu : slack variables for inequality constraints
            alpha : continuation parameter
        '''

        o = self
        
        # y = [x, lambda, gamma]
        nx = len(o.x) # number of state variables
        np = len(o.p) # number of paraters
        y = [0] * (2*nx + np)
        for i in range(nx):
            y[i] = o.x[i]
            y[i+nx] = Symbol('_lambda%d' % (i+1))
        for i in range(np):
            y[2*nx+i] = Symbol('_gamma%d' % (i+1))
        
        # z = [u, theta, mu, nu]
        nu = len(o.u) # number of control varaibles
        nc = len(o.c) # number of equality constraints
        nd = len(o.d) # number of inequality constraints
        z = [0] * (nu + nc + 2*nd)
        for i in range(nu):
            z[i] = o.u[i]
        for i in range(nc):
            z[nu+i] = Symbol('_theta%d' % (i+1))
        for i in range(nd):
            z[nu+nc+i] = Symbol('_mu%d' % (i+1))
            z[nu+nc+nd+i] = Symbol('_nu%d' % (i+1))
        
        # p = [p, nu_i, nu_f]
        ngamma = len(o.Gamma) # number of initial constraints
        npsi = len(o.Psi) # number of final constraints
        p = [0] * (np + ngamma + npsi)
        for i in range(np):
            p[i] = o.p[i]
        for i in range(ngamma):
            p[np+i] = Symbol('_kappa_i%d' % (i+1))
        for i in range(npsi):
            p[np+ngamma+i] = Symbol('_kappa_f%d' % (i+1))

        # alpha : continuation parameter
        _alpha = Symbol('_alpha')

        # form the Hamiltonian
        H = self.L
        for i in range(len(o.x)):
            H = H + y[nx+i]*o.f[i]
        for i in range(len(o.c)):
            H = H + z[nu+i]*o.c[i]
        for i in range(len(o.d)):
            H = H + z[nu+nc+i]*o.d[i]
        
        # print H
        
        # form the states, co-states and parameters : f
        # 
        # f = [H_lambda; -H_x; H_p]
        #
        f = [0] * (2*nx + np)
        for i in range(nx):
            f[i] = diff(H, y[nx+i]) # = o.f[i]
            f[nx+i] = -diff(H, o.x[i])
        for i in range(np):
            f[2*nx+i] = diff(H, o.p[i])
        
        # form the stationary condition and complementarity conditions : g
        #
        # g = [H_u; H_theta; H_mu+nu; mu*nu]
        #
        g = [0] * (nu + nc + 2*nd)
        for i in range(nu):
            g[i] = diff(H, o.u[i])
        for i in range(nc):
            g[nu+i] = diff(H, z[nu+i]) # =  o.c[i]
        for i in range(nd):
            g[nu+nc+i] = diff(H, z[nu+nc+i]) + z[nu+nc+nd+i] # = o.d[i] + z[nu+nc+nd+i]
            # g[nu+nc+nd+i] = z[nu+nc+i]*z[nu+nc+nd+i]
            g[nu+nc+nd+i] = z[nu+nc+i] + z[nu+nc+nd+i] - sqrt(z[nu+nc+i]**2 + z[nu+nc+nd+i]**2 + 2 * _alpha) 
        
        # form the boundary conditions : r
        #
        # r_initial = [Gamma; lambda + Gamma_x*nu_i; gamma]
        # f_final = [Psi; lambda - phi_x - Psi_x*nu_f; gamma + phi_p + Psi_p*nu_f]
        #
        ri = [0] * (ngamma + nx + np)
        for i in range(ngamma):
            ri[i] = o.Gamma[i]
        Gamma_x_nu_i = [0] * nx
        if ngamma > 0:
            for i in range(nx):
                for j in range(ngamma):
                    Gamma_x_nu_i[i] += diff(o.Gamma[j], y[i]) * p[np+j]
        for i in range(nx):
            ri[ngamma+i] = y[nx+i] + Gamma_x_nu_i[i]
        for i in range(np):
            ri[ngamma+nx+i] = y[2*nx+i]
        
        rf = [0] * (npsi + nx + np)
        for i in range(npsi):
            rf[i] = o.Psi[i]
        Psi_x_nu_f = [0] * nx
        if npsi > 0:
            for i in range(nx):
                for j in range(npsi):
                    Psi_x_nu_f[i] += diff(o.Psi[j], y[i]) * p[np+ngamma+j]
        for i in range(nx):
            rf[npsi+i] = y[nx+i] - diff(o.phi, y[i]) - Psi_x_nu_f[i]
        if np > 0:
            Psi_p_nu_f = [0] * np
            for i in range(np):
                for j in range(npsi):
                    Psi_p_nu_f[i] += diff(o.Psi[j], p[i]) * p[np+ngamma+j]
            for i in range(np):
                rf[npsi+nx+i] = y[2*nx+i] + diff(o.phi, p[i]) + Psi_p_nu_f[i]

        o.write_abvp_c_code(y, z, p, f, g, ri, rf, len(o.d))

    def write_abvp_c_code(self, y, z, p, f, g, ri, rf, n_ineq):
        o = self
        '''
            construct a BVP_DAE class
        '''
        
        # constant values defined in the OCP file
        constants = ''
        if o.constants:
            for i in range(len(o.constants)):
                constants += '\t\t{0} {1}\n'.format(o.constants[i], o.constants_value[i])
        
        # f(y,z,p)
        fg_vars = ''
        for i in range(len(y)):
            fg_vars += '\t\t{0} = _y[{1}]\n'.format(y[i], i)
        for i in range(len(z)):
            fg_vars += '\t\t{0} = _z[{1}]\n'.format(z[i], i)
        for i in range(len(o.p)):
            fg_vars += '\t\t{0} = _p[{1}]\n'.format(o.p[i], i)
            
        abvp_f = '\tdef _abvp_f(self, _y, _z, _p, _f):\n'
        if o.constants:
            abvp_f += constants
        abvp_f += fg_vars
        for i in range(len(f)):
            abvp_f += '\t\t_f[{0}] = {1}\n'.format(i, ccode(f[i]))
        abvp_f += '\n'

        # g(y,z,p)
        abvp_g = '\tdef _abvp_g(self, _y, _z, _p, _alpha, _g):\n'
        if o.constants:
            abvp_g += constants
        abvp_g += fg_vars
        for i in range(len(g)):
            abvp_g += '\t\t_g[{0}] = {1}\n'.format(i, ccode(g[i]))
        abvp_g += '\n'
        
        # r(y0, y1, p)
        ri_vars = ''
        for i in range(len(p)):
            ri_vars += '\t\t{0} = _p[{1}]\n'.format(p[i], i)
        for i in range(len(y)):
            ri_vars += '\t\t{0} = _y0[{1}]\n'.format(y[i], i)
        abvp_r = '\tdef _abvp_r(self, _y0, _y1, _p, _r):\n'
        if o.constants:
            abvp_r += constants
        abvp_r += ri_vars
        abvp_r += '\t\t# initial conditions\n'
        for i in range(len(ri)):
            abvp_r += '\t\t_r[{0}] = {1}\n'.format(i, ccode(ri[i]))
        
        rf_vars = ''
        for i in range(len(y)):
            rf_vars += '\t\t{0} = _y1[{1}]\n'.format(y[i], i)
        nri = len(ri)
        abvp_r += '\t\t# final conditions\n'
        abvp_r += rf_vars
        for i in range(len(rf)):
            abvp_r += '\t\t_r[{0}] = {1}\n'.format(i+nri, ccode(rf[i]))
        abvp_r += '\n'
        
        # df(y, z, p)
        abvp_Df = '\tdef _abvp_Df(self, _y, _z, _p, _Df):\n'
        if o.constants:
            abvp_Df += constants
        abvp_Df += fg_vars
        ny = len(y)
        nz = len(z)
        np = len(p)
        for i in range(len(f)):
            for j in range(len(y)):
                df = diff(f[i], y[j])
                if df != 0:
                    abvp_Df += '\t\t_Df[{0}][{1}] = {2}\n'.format(i, j, ccode(df))
            for j in range(len(z)):
                df = diff(f[i], z[j])
                if df != 0:
                    abvp_Df += '\t\t_Df[{0}][{1}] = {2}\n'.format(i, ny+j, ccode(df))
            for j in range(len(p)):
                df = diff(f[i], p[j])
                if df != 0:
                    abvp_Df += '\t\t_Df[{0}][{1}] = {2}\n'.format(i, ny+nz+j, ccode(df))
        abvp_Df += '\n'
        
        # dg(y, z, p)
        abvp_Dg = '\tdef _abvp_Dg(self, _y, _z, _p, _alpha, _Dg):\n'
        if o.constants:
            abvp_Dg += constants
        abvp_Dg += fg_vars
        for i in range(len(g)):
            for j in range(len(y)):
                df = diff(g[i], y[j])
                if df != 0:
                    abvp_Dg += '\t\t_Dg[{0}][{1}] = {2}\n'.format(i, j, ccode(df))
            for j in range(len(z)):
                df = diff(g[i], z[j])
                if df != 0:
                    abvp_Dg += '\t\t_Dg[{0}][{1}] = {2}\n'.format(i, ny+j, ccode(df))
            for j in range(len(p)):
                df = diff(g[i], p[j])
                if df != 0:
                    abvp_Dg += '\t\t_Dg[{0}][{1}] = {2}\n'.format(i, ny+nz+j, ccode(df))
        abvp_Dg += '\n'
        
        # dr(y0, y1, p)
        abvp_Dr = '\tdef _abvp_Dr(self, _y0, _y1, _p, _Dr):\n'
        if o.constants:
            abvp_Dr += constants
        abvp_Dr += ri_vars
        abvp_Dr += '\t\t# initial conditions\n'
        for i in range(len(ri)):
            for j in range(len(y)):
                dr = diff(ri[i], y[j])
                if dr != 0:
                    abvp_Dr += '\t\t_Dr[{0}][{1}] = {2}\n'.format(i, j, ccode(dr))
            for j in range(len(p)):
                dr = diff(ri[i], p[j])
                if dr != 0:
                    abvp_Dr += '\t\t_Dr[{0}][{1}] = {2}\n'.format(i, 2*ny+j, ccode(dr))
        nri = len(ri)
        abvp_Dr += '\t\t# final conditions\n'
        abvp_Dr += rf_vars
        for i in range(len(rf)):
            for j in range(len(y)):
                dr = diff(rf[i], y[j])
                if dr != 0:
                    abvp_Dr += '\t\t_Dr[{0}][{1}] = {2}\n'.format(i+nri, ny+j, ccode(dr))
            for j in range(len(p)):
                dr = diff(rf[i], p[j])
                if dr != 0:
                    abvp_Dr += '\t\t_Dr[{0}][{1}] = {2}\n'.format(i+nri, 2*ny+j, ccode(dr))
        abvp_Dr += '\n'
        
        # main()
        abvp_main = '\tdef __init__(self) :\n'
        abvp_main += '\t\tself.size_y = {0}\n'.format(len(y))
        abvp_main += '\t\tself.size_z = {0}\n'.format(len(z))
        abvp_main += '\t\tself.size_p = {0}\n'.format(len(p))
        abvp_main += '\t\tself.size_inequality = {0}\n'.format(n_ineq)
        '''
        abvp_main += '\tABVPDAE _bvp = ABVPDAENew(_ny, _nz, _np, _ninequality);\n'
        abvp_main += '\t_bvp->f = _abvp_f;\n'
        abvp_main += '\t_bvp->g = _abvp_g;\n'
        abvp_main += '\t_bvp->r = _abvp_r;\n'
        abvp_main += '\t_bvp->Df = _abvp_Df;\n'
        abvp_main += '\t_bvp->Dg = _abvp_Dg;\n'
        abvp_main += '\t_bvp->Dr = _abvp_Dr;\n'
        abvp_main += '\tMSHOOTDAE _m = MSHOOTDAENew();\n'
        abvp_main += '\t_m->bvp = _bvp;\n'
        '''
        abvp_main += '\t\tself.tolerance = {0}\n'.format(o.tolerance)
        abvp_main += '\t\tself.maximum_nodes = {0}\n'.format(o.maximum_nodes)
        abvp_main += '\t\tself.maximum_newton_iterations = {0}\n'.format(o.maximum_newton_iterations)
        abvp_main += '\t\tself.maximum_mesh_refinements = {0}\n'.format(o.maximum_mesh_refinements)
        # abvp_main += '\t_m->display = {0};\n'.format(o.display)
        if o.input_file != "":
            abvp_main += '\tMSHOOTDAEData _data = MSHOOTDAEReadData("{0}");\n'.format(o.input_file)
            abvp_main += '\tif (_data.error != 0) { RuntimeWarning("Unable to read input file"); ABVPDAEDelete(_bvp); MSHOOTDAEDelete(_m); return _data.error; }\n'
            abvp_main += '\tVector _T0 = _data.T;\n'
            #abvp_main += '\tMatrix _Y0 = _data.Y;\n'
            #abvp_main += '\tMatrix _Z0 = _data.Z;\n'
            #abvp_main += '\tVector _P0 = NULL; if (_data.P != NULL) { _P0 = _data.P; }\n'
            abvp_main += '\tMatrix _Y0 = NULL; if (_ny > 0) { _Y0 = MatrixNew(_T0->r, _ny); MatrixSetAllTo(_Y0, 1.0); }\n'
            abvp_main += '\tMatrix _Z0 = NULL; if (_nz > 0) { _Z0 = MatrixNew(_T0->r, _nz); MatrixSetAllTo(_Z0, 10.0); }\n'
            abvp_main += '\tVector _P0 = NULL; if (_np > 0) { _P0 = VectorNew(_np); VectorSetAllTo(_P0, 10.0); }\n'
            abvp_main += '\t_pack_YZP(_Y0, _Z0, _P0, _data);\n'
            abvp_main += '\tif (_data.Y != NULL) MatrixDelete(_data.Y);\n'
            abvp_main += '\tif (_data.Z != NULL) MatrixDelete(_data.Z);\n'
            abvp_main += '\tif (_data.P != NULL) VectorDelete(_data.P);\n'
            if o.state_estimate or o.control_estimate or o.parameter_estimate:
                abvp_main += '\t_solution_estimate(_T0, _Y0, _Z0, _P0);\n'
        else:
            abvp_main += '\t\tself.N = {0}\n'.format(o.nodes)
            abvp_main += '\t\tself.t_initial = {0}\n'.format(o.t_i)
            abvp_main += '\t\tself.t_final = {0}\n'.format(o.t_f)
            abvp_main += '\t\tself.T0 = np.linspace(self.t_initial, self.t_final, self.N)\n'
            abvp_main += '\t\tself.Y0 = np.ones((self.N, self.size_y), dtype = np.float64)\n'
            abvp_main += '\t\tself.Z0 = np.ones((self.N, self.size_z), dtype = np.float64)\n'
            if np > 0:
                abvp_main += '\t\tself.P0 = np.ones((self.size_p), dtype = np.float64)\n'
            else:
                abvp_main += '\t\tself.P0 = np.ones((0), dtype = np.float64)\n'
            '''
            abvp_main += '\tMatrixSetAllTo(_Y0, 1.0)\n'
            abvp_main += '\tMatrixSetAllTo(_Z0, 1.0)\n'
            if np > 0:
                abvp_main += '\tVectorSetAllTo(_P0, 1.0)\n'
            '''
            if o.state_estimate or o.control_estimate or o.parameter_estimate:
                abvp_main += '\t\tself._solution_estimate(self.T0, self.Y0, self.Z0, self.P0)\n'

        # initial solution estimate
        if o.state_estimate or o.control_estimate or o.parameter_estimate:
            abvp_main += '\n'
            abvp_main += '\tdef _solution_estimate(self, _T, _Y, _Z, _P):\n'
            # abvp_header += '\tint i;\n'
            # abvp_header += '\tdouble t;\n'
            abvp_main += '\t\tN = _T.shape[0]\n'
            abvp_main += '\t\tfor i in range(N):\n'
            abvp_main += '\t\t\tt = _T[i];\n'
            if o.state_estimate:
                for j in range(len(o.state_estimate)):
                    abvp_main += '\t\t\t_Y[i][{0}] = {1};\n'.format(j, o.state_estimate[j])
            if o.control_estimate:
                for j in range(len(o.control_estimate)):
                    #abvp_header += '\t\t_Z->e[i][{0}] = {1};\n'.format(j, ccode(o.control_estimate[j]))
                    abvp_main += '\t\t\t_Z[i][{0}] = {1};\n'.format(j, o.control_estimate[j])
            #abvp_header += '\t\t}\n'
            abvp_main += '\n'
            abvp_main += '\t\tif (_P.shape[0] != 0):\n'
            abvp_main += '\t\t\tfor i in range(self.size_p):\n'
            if o.parameter_estimate:
                for j in range(len(o.parameter_estimate)):
                    abvp_main += '\t\t\t\t_P[{0}] = {1};\n'.format(j, o.parameter_estimate[j])
            else:
                abvp_main += '\t\t\t\t_P0 = np.ones((self.size_p), dtype = np.float64)\n'
        
        '''
        abvp_main += '\tint _err = MSHOOTDAESolve(_m, _T0, _Y0, _Z0, _P0);\n'
        abvp_main += '\tprintf("MSHOOTDAESolve exit code: %d\\n", _err);\n'
        abvp_main += '\t_err = MSHOOTDAEWriteData(_m, "{0}");\n'.format(o.output_file)
        abvp_main += '\tABVPDAEDelete(_bvp);\n\tMSHOOTDAEDelete(_m);\n\tVectorDelete(_T0);\n\tMatrixDelete(_Y0);\n\tMatrixDelete(_Z0);\n'
        if np > 0:
            abvp_main += '\tVectorDelete(_P0);\n'
        abvp_main += '\treturn _err;\n}\n'
        '''
        
        # header, import necessary packages and class header
        abvp_header = '# Created by OCP.py\n'
        abvp_header += 'from math import *\n'
        abvp_header += 'import numpy as np\n\n'
        abvp_header += 'class bvp_dae:\n'

        '''
        abvp_header += '#include <stdio.h>\n'
        abvp_header += '#include <float.h>\n'
        abvp_header += '#include <math.h>\n'
        abvp_header += '#include <stdlib.h>\n'
        abvp_header += '#include "runtime.h"\n'
        abvp_header += '#include "vector.h"\n'
        abvp_header += '#include "matrix.h"\n'
        abvp_header += '#include "lu.h"\n'
        abvp_header += '#include "qr.h"\n'
        abvp_header += '#include "abd.h"\n'
        abvp_header += '#include "equations.h"\n'
        abvp_header += '#include "newtonsmethod.h"\n'
        abvp_header += '#include "abvpdae.h"\n'
        abvp_header += '#include "mshootdae.h"\n\n'

        abvp_header += 'double _rho_ = 1.0e3;\n\n'
        abvp_header += 'double _mu_ = 1.0e-1;\n\n'
        if o.constants:
            for i in range(len(o.constants)):
                abvp_header += 'double {0} {1};\n'.format(o.constants[i], o.constants_value[i])
        
        # pack data from the input file
        if o.input_file != "":
            abvp_header += 'void _pack_YZP(Matrix _Y, Matrix _Z, Vector _P, MSHOOTDAEData _data) {\n'
            abvp_header += '    int _i, _j, _n, _m;\n'
            abvp_header += '    if ((_Y != NULL) && (_data.Y != NULL)) {\n'
            abvp_header += '        _n = _data.Y->r < _Y->r ? _data.Y->r : _Y->r;\n'
            abvp_header += '        _m = _data.Y->c < _Y->c ? _data.Y->c : _Y->c;\n'
            abvp_header += '        for (_i = 0; _i < _n; _i++) {\n'
            abvp_header += '            for (_j = 0; _j < _m; _j++) {\n'
            abvp_header += '                _Y->e[_i][_j] = _data.Y->e[_i][_j];\n'
            abvp_header += '            }\n'
            abvp_header += '        }\n'
            abvp_header += '    }\n'
            abvp_header += '    if ((_Z != NULL) && (_data.Z != NULL)) {\n'
            abvp_header += '        _n = _data.Z->r < _Z->r ? _data.Z->r : _Z->r;\n'
            #abvp_header += '        _m = _data.Z->c < _Z->c ? _data.Z->c : _Z->c;\n'
            abvp_header += '        _m = _data.Z->c < {0} ? _data.Z->c : {1}; // only read in enough dtat to fill the controls\n'.format(len(o.u), len(o.u))
            #abvp_header += '        _m = {0};\n'.format(len(o.u))
            abvp_header += '        for (_i = 0; _i < _n; _i++) {\n'
            abvp_header += '            for (_j = 0; _j < _m; _j++) {\n'
            abvp_header += '                _Z->e[_i][_j] = _data.Z->e[_i][_j];\n'
            abvp_header += '            }\n'
            abvp_header += '        }\n'
            abvp_header += '    }\n'
            abvp_header += '    if ((_P != NULL) && (_data.P != NULL)) {\n'
            abvp_header += '        _n = _data.P->r < _P->r ? _data.P->r : _P->r;\n'
            abvp_header += '        for (_i = 0; _i < _n; _i++) {\n'
            abvp_header += '            _P->e[_i] = _data.P->e[_i];\n'
            abvp_header += '        }\n'
            abvp_header += '    }\n'
            abvp_header += '}\n'
        # initial solution estimate
        if o.state_estimate or o.control_estimate or o.parameter_estimate:
            abvp_header += 'void _solution_estimate(Vector _T, Matrix _Y, Matrix _Z, Vector _P) {\n'
            abvp_header += '\tint i;\n'
            abvp_header += '\tdouble t;\n'
            abvp_header += '\tfor (i = 0; i < _T->r; i++) {\n'
            abvp_header += '\t\tt = _T->e[i];\n'
            if o.state_estimate:
                for j in range(len(o.state_estimate)):
                    abvp_header += '\t\t_Y->e[i][{0}] = {1};\n'.format(j, o.state_estimate[j])
            if o.control_estimate:
                for j in range(len(o.control_estimate)):
                    #abvp_header += '\t\t_Z->e[i][{0}] = {1};\n'.format(j, ccode(o.control_estimate[j]))
                    abvp_header += '\t\t_Z->e[i][{0}] = {1};\n'.format(j, o.control_estimate[j])
            #abvp_header += '\t\t}\n'
            abvp_header += '\t}\n'
            abvp_header += '\tif (_P != NULL) {\n'
            abvp_header += '\t\tfor (i = 0; i < _P->r; i++) {\n'
            if o.parameter_estimate:
                for j in range(len(o.parameter_estimate)):
                    abvp_header += '\t\t\t_P->e[{0}] = {1};\n'.format(j, o.parameter_estimate[j])
            abvp_header += '\t\t}\n'
            abvp_header += '\t}\n'
            abvp_header += '}\n'
        '''
        print (abvp_header)
        print (abvp_main)
        print (abvp_f)
        print (abvp_g)
        print (abvp_r)
        print (abvp_Df)
        print (abvp_Dg)
        print (abvp_Dr)
        # print (abvp_main)

def _make_variables_py0(line, typ):
    rline = ''
    k = 0
    np = line.count(',') + 1
    rline += '_o.%s = [0] * %d\n' % (typ,np)
    i = line.find('[')
    while (line[i] != ']'):
        p = ''
        while ((line[i + 1] != ',') \
            and (line[i + 1] != ']')):
            p += line[i + 1]
            i += 1
        rline += '%s = Symbol(\'%s\')\n' % (p, p)
        rline += '_o.%s[%d] = %s\n' % (typ, k, p)
        i += 1
        k += 1

    return (rline, k)
    
def _make_constants_py0(line):
    global _raw_line_number_
    rline = ''
    k = 0
    np = line.count(',') + 1
    rline += '_o.constants = [0] * %d\n' % np
    rline += '_o.constants_value = [0] * %d\n' % np
    i = line.find('[')
    try:
        while (line[i] != ']'):
            p = ''
            pv = ''
            while (line[i + 1] != '='): # get the constants
                p += line[i + 1]
                i += 1
            while ((line[i + 1] != ',')
                and (line[i + 1] != ']')): # get the value
                pv += line[i + 1]
                i += 1
            rline += '%s = Symbol(\'%s\')\n' % (p, p)
            rline += '_o.constants[%d] = \'%s\'\n' % (k, p)
            rline += '_o.constants_value[%d] = \'%s\'\n' % (k, pv)
            i += 1
            k += 1
    except:
        print ('Syntax error on line: ', _raw_line_number_)
        print ('Exception: ', sys.exc_type, sys.exc_value)
        raise SyntaxError

    return rline

def _line_to_ocp(line):
    global _n_, _m1_, _m2_, _m3_, _m4_
    rline = ''
    if (line.find('Constants') == 0):
        rline = '# ' + line + '\n'
        vline = _make_constants_py0(line)
        rline += vline
    elif (line.find('StateVariables') == 0):
        rline = '# ' + line + '\n'
        (vline, _n_) = _make_variables_py0(line, 'x')
        rline += vline
    elif (line.find('ControlVariables') == 0):
        rline = '# ' + line + '\n'
        (vline, _n_) = _make_variables_py0(line, 'u')
        rline += vline
    elif (line.find('ParameterVariables') == 0):
        rline = '# ' + line + '\n'
        (vline, _n_) = _make_variables_py0(line, 'p')
        rline += vline
    elif (line.find('TerminalPenalty') == 0):
        rline = line.replace('TerminalPenalty', '_o.phi')
    elif (line.find('CostFunctional') == 0):
        rline = line.replace('CostFunctional', '_o.L')
    elif (line.find('InitialConstraints') == 0):
        rline = line.replace('InitialConstraints', '_o.Gamma')
    elif (line.find('TerminalConstraints') == 0):
        rline = line.replace('TerminalConstraints', '_o.Psi')
    elif (line.find('DifferentialEquations') == 0):
        rline = line.replace('DifferentialEquations', '_o.f')
    elif (line.find('InequalityConstraints') == 0):
        rline = line.replace('InequalityConstraints', '_o.d')
    elif (line.find('EqualityConstraints') == 0):
        rline = line.replace('EqualityConstraints', '_o.c')
    elif (line.find('InitialTime') == 0):
        rline = line.replace('InitialTime', '_o.t_i')
    elif (line.find('FinalTime') == 0):
        rline = line.replace('FinalTime', '_o.t_f')
    elif (line.find('Nodes') == 0):
        rline = line.replace('Nodes', '_o.nodes')
    elif (line.find('Tolerance') == 0):
        rline = line.replace('Tolerance', '_o.tolerance')
    elif (line.find('InputFile') == 0):
        rline = line.replace('InputFile', '_o.input_file')
    elif (line.find('OutputFile') == 0):
        rline = line.replace('OutputFile', '_o.output_file')
    elif (line.find('Display') == 0):
        rline = line.replace('Display', '_o.display')
    elif (line.find('MaximumMeshRefinements') == 0):
        rline = line.replace('MaximumMeshRefinements', '_o.maximum_mesh_refinements')
    elif (line.find('MaximumNewtonIterations') == 0):
        rline = line.replace('MaximumNewtonIterations', '_o.maximum_newton_iterations')
    elif (line.find('MaximumNodes') == 0):
        rline = line.replace('MaximumNodes', '_o.maximum_nodes')
    elif (line.find('StateEstimate') == 0):
        rline = line.replace('StateEstimate', '_o.state_estimate')
    elif (line.find('ControlEstimate') == 0):
        rline = line.replace('ControlEstimate', '_o.control_estimate')
    elif (line.find('ParameterEstimate') == 0):
        rline = line.replace('ParameterEstimate', '_o.parameter_estimate')
    else:
        rline = line
    return rline

def _ocp_translate(inpt, typ):
    """
        Translate the OCP script file
        If typ == 0 -> C MSHOOTDAE
        If typ == 1 -> C RKOCP
        If typ == 2 -> C ROWOCP
        If typ == 3 -> FORTRAN TOMP
    """
    global _raw_line_number_
    fid = open(inpt, 'r')
    s = ''
    rawline = ''
    # raw lines from the input file
    _raw_lines_ = ''
    while 1:
        line = fid.readline()
        _raw_lines_ += line
        # print (_raw_lines_)
        _raw_line_number_ += 1
        if not line:
            break
        indx = line.find('#')
        if (indx < 0):
            rline = line
        else:
            rline = line[0 : indx]
        for i in range(len(rline)):
            # ignore newline, return and white spaces
            if ((rline[i] != '\n') \
                and (rline[i] != '\r') \
                and (rline[i] != '\t') \
                and (rline[i] != ' ')):
                rawline += rline[i]
                if (rline[i] == ';'):
                    #print rawline
                    rawline = _line_to_ocp(rawline)
                    #print rawline
                    s += rawline + '\n'
                    rawline = ''
    r = ''
    r += 'from OCP2ABVP import *\n'
    r += 'NO = 0\n'
    r += 'YES = 1\n'
    r += 't = Symbol(\'t\')\n'
    r += '_o = OCP()\n'
    s += '_o.make_abvp()\n'
    t = '%s\n%s' % (r, s)
    '''
    print('r:')
    print(r)
    print('r done.')
    print('s:')
    print(s)
    print('s done.')
    print('t')
    print(t)
    print('t done')
    print('Raw lines:')
    print(_raw_lines_)
    print('raw lines done!')
    '''
    exec(t)
    print ('\'\'\'\n')
    print (_raw_lines_)
    print ('\'\'\'\n')

if __name__ == '__main__':
    """
    u1 = Symbol('u1')
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    o = OCP()
    o.x = [x1, x2]
    o.u = [u1]
    o.L = 0.5*u1*u1
    o.f = [x2, u1]
    o.Gamma = [x1-1, x2-1]
    o.Psi = [x1, x2]
    """
    
    """
    u1 = Symbol('u1')
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    rho = Symbol('rho')
    o = OCP()
    o.x = [x1, x2]
    o.u = [u1]
    o.L = -x2 + 0.5*u1*u1*rho
    o.f = [x2, u1]
    o.d = [u1 - 1.0, -1.0 - u1]
    o.Gamma = [x1-1.0, x2-1.0]
    o.Psi = [x2]
    #o.input_file = "ex3.data"
    o.output_file = "ex4.data"
    """
    """
    u1 = Symbol('u1')
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    o = OCP()
    o.x = [x1, x2]
    o.u = [u1]
    o.Gamma = [x1-0.05, x2]
    o.phi = x1*x1 + x2*x2
    o.L = 0.5*(x1*x1 + x2*x2 + 0.1*u1*u1)
    a1 = x1 + 0.25
    a2 = x2 + 0.5
    a3 = x1 + 2.0
    a4 = a2*exp(25.0*x1/a3)
    o.f = [-2.0*a1 + a4 - a1*u1, 0.5 - x2 - a4]
    o.output_file = "ex7a.data"
    o.t_i = 0
    o.t_f = 0.78
    o.nodes = 79
    
    o.make_abvp()
    """
    try:
        ocp = sys.argv[1]
    except:
        print ('Unable to read input file')
        print ('Exception: ', sys.exc_type, sys.exc_value)

    _ocp_translate(ocp, 0)
