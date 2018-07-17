
import numpy as np

class lobatto:
    '''
        Define the lobatto coefficients for collocation points
    '''

    '''
        Input: m - stages of the lobatto coefficients
    '''
    def __init__(self, m):
        self.stages = m
        self.order = m + 1
        self.A = np.zeros((m, m), dtype = np.float64)
        self.b = np.zeros((m), dtype = np.float64)
        self.c = np.zeros((m), dtype = np.float64)
        self.d = np.zeros((m + 1), dtype = np.float64)
        self.w = np.zeros((m + 1), dtype = np.float64)

        if m == 3:
            self.A[0, 0] = 0.0000000000000000e+00;
            self.A[0, 1] = 0.0000000000000000e+00;
            self.A[0, 2] = 0.0000000000000000e+00;
            self.A[1, 0] = 2.0833333333333337e-01;
            self.A[1, 1] = 3.3333333333333337e-01;
            self.A[1, 2] = -4.1666666666666671e-02;
            self.A[2, 0] = 1.6666666666666652e-01;
            self.A[2, 1] = 6.6666666666666674e-01;
            self.A[2, 2] = 1.6666666666666663e-01;
            self.b[0] = 1.6666666666666652e-01;
            self.b[1] = 6.6666666666666674e-01;
            self.b[2] = 1.6666666666666663e-01;
            self.c[0] = 0.0000000000000000e+00;
            self.c[1] = 5.0000000000000000e-01;
            self.c[2] = 1.0000000000000000e+00;
            self.w[0] = 8.3333333333333925e-02;
            self.w[1] = 4.1666666666666519e-01;
            self.w[2] = 4.1666666666666607e-01;
            self.w[3] = 8.3333333333333259e-02;
            self.d[0] = 0.0000000000000000e+00;
            self.d[1] = 2.7639320225002106e-01;
            self.d[2] = 7.2360679774997894e-01;
            self.d[3] = 1.0000000000000000e+00;
