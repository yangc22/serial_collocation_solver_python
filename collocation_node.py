import numpy as np

class collocation_node:
    """
        Node class for collocation solver, save the required DAE variables
        at each time node
    """

    '''
        Input: size_y - size of the ODE variables
               size_z - size of the DAE variables
               size_p - size of the parameters
               m: stages of the lobatto coefficients 
        Construct the empty collocation node with all the 
        fields set to zero.
    '''
    def __init__(self, size_y, size_z, size_p, m):
        self.size_y = size_y
        self.size_z = size_z
        self.size_p = size_p
        self.m = m
        self.delta_t = []
        self.tspan = []
        self.y = np.zeros((size_y), dtype = np.float64)
        self.z = np.zeros((size_z), dtype = np.float64)
        self.y_dot = np.zeros((size_y, m), dtype = np.float64)
        self.z_tilda = np.zeros((size_z, m), dtype = np.float64)
        self.p = np.zeros((size_p), dtype = np.float64)
        self.y_tilda = np.zeros((size_y, m), dtype = np.float64)
        self.f_a = np.zeros(((size_y + size_z) * m), dtype = np.float64)
        self.f_b = np.zeros((size_y), dtype = np.float64)
        self.J = np.zeros(((size_y + size_z) * m, size_y), dtype = np.float64)
        self.W = np.zeros(((size_y + size_z) * m, (size_y + size_z) * m), dtype = np.float64)
        self.V = np.zeros(((size_y + size_z) * m, size_p), dtype = np.float64)
        self.D = np.zeros((size_y, (size_y + size_z) * m), dtype = np.float64)
        self.A = np.zeros((size_y, size_y), dtype = np.float64)
        self.C = np.zeros((size_y, size_y), dtype = np.float64)
        self.H = np.zeros((size_y, size_p), dtype = np.float64)
        self.b = np.zeros((size_y), dtype = np.float64)
        '''
        self.C_Tilde = np.zeros()
        self.G_Tilde = np.zeros()
        self.H_Tilde = np.zeros()
        self.b_Tilde = np.zeros()
        self.R = np.zeros()
        self.E = np.zeros()
        self.G = np.zeros()
        self.K = np.zeros()
        self.d = np.zeros()
        self.Rp = np.zeros()
        self.dp = np.zeros()
        self.delta_p = np.zeros()
        self.delta_y = np.zeros()
        self.delta_k = np.zeros()
        '''

    '''
        Input: y - a size_y vector of the value of the ODE variables
        Set the ODE variable y.
    '''
    def set_y(self, y):
        for i in range(self.size_y):
            self.y[i] = y[i]

    '''
        Input: z - a size_z vector of the value of the DAE variables
        Set the DAE variable z.
    '''
    def set_z(self, z):
        for i in range(self.size_z):
            self.z[i] = z[i]

    '''
        Input: p - a size_p x 1 vector of the value of the parameter variables
        Set the parameter variable p.
    '''
    def set_p(self, p):
        for i in range(self.size_p):
            self.p[i] = p[i]

    '''
        Input: delta_t - a double representing the interval of the time 
                         span between the current node and the next node
        Set the time interval delta_t.
    '''
    def set_delta_t(self, delta_t):
        self.delta_t = delta_t

    '''
        Input: tspan - a 1 x 2 vector representing the time span between the current node and the next node
        Set the time interval tspan.
    '''
    def set_tspan(self, tspan):
        self.tspan = tspan

    def set_y_dot(self, y_dot, j):
        for i in range(self.size_y):
            self.y_dot[i][j] = y_dot[i]

    def set_z_tilda(self, z_tilda, j):
        for i in range(self.size_z):
            self.z_tilda[i][j] = z_tilda[i]

