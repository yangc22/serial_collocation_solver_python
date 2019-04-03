import numpy as np

# Define the gaussian coefficients for collocation points


class gauss:

    # Input: n - stages of the gaussian coefficients
    def __init__(self, n):
        self.n = n

        if n == 2:
            x1 = 0.577350269189626
            w1 = 1.000000000000000
            t = [(-x1+1)*0.5, (x1+1)*0.5]
            w = [w1*0.5, w1*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 3:
            x0 = 0.000000000000000
            x1 = 0.774596669241483
            w0 = 0.888888888888889
            w1 = 0.555555555555556
            t = [(-x1+1)*0.5, 0.5, (x1+1)*0.5]
            w = [w1*0.5, w0*0.5, w1*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 4:
            x1 = 0.339981043584856
            x2 = 0.861136311594053
            w1 = 0.652145154862546
            w2 = 0.347854845137454
            t = [(-x2+1)*0.5, (-x1+1)*0.5, (x1+1)*0.5, (x2+1)*0.5]
            w = [w2*0.5, w1*0.5,  w1*0.5, w2*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 5:
            x0 = 0.000000000000000
            x1 = 0.538469310105683
            x2 = 0.906179845938664
            w0 = 0.568888888888889
            w1 = 0.478628670499366
            w2 = 0.236926885056189
            t = [(-x2+1)*0.5, (-x1+1)*0.5, 0.5, (x1+1)*0.5, (x2+1)*0.5]
            w = [w2*0.5, w1*0.5, w0*0.5, w1*0.5, w2*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 6:
            x1 = 0.238619186083197
            x2 = 0.661209386466265
            x3 = 0.932469514203152
            w1 = 0.467913934572691
            w2 = 0.360761573048139
            w3 = 0.171324492379170
            t = [(-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5]
            w = [w3*0.5, w2*0.5, w1*0.5,  w1*0.5, w2*0.5, w3*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 7:
            x0 = 0.000000000000000
            x1 = 0.405845151377397
            x2 = 0.741531185599394
            x3 = 0.949107912342759
            w0 = 0.417959183673469
            w1 = 0.381830050505119
            w2 = 0.279705391489277
            w3 = 0.129484966168870
            t = [(-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, 0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5]
            w = [w3*0.5, w2*0.5, w1*0.5,  w0*0.5, w1*0.5, w2*0.5, w3*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 8:
            x1 = 0.183434642495650
            x2 = 0.525532409916329
            x3 = 0.796666477413627
            x4 = 0.960289856497536
            w1 = 0.362683783378362
            w2 = 0.313706645877887
            w3 = 0.222381034453374
            w4 = 0.101228536290376
            t = [(-x4+1)*0.5, (-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5, (x4+1)*0.5]
            w = [w4*0.5, w3*0.5, w2*0.5, w1*0.5,  w1*0.5, w2*0.5, w3*0.5, w4*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 9:
            x0 = 0.000000000000000
            x1 = 0.324253423403809
            x2 = 0.613371432700590
            x3 = 0.836031107326636
            x4 = 0.968160239507626
            w0 = 0.330239355001260
            w1 = 0.312347077040003
            w2 = 0.260610696402935
            w3 = 0.180648160694857
            w4 = 0.081274388361574
            t = [(-x4+1)*0.5, (-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, 0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5, (x4+1)*0.5]
            w = [w4*0.5, w3*0.5, w2*0.5, w1*0.5,  w0*0.5, w1*0.5, w2*0.5, w3*0.5, w4*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 10:
            x1 = 0.148874338981631
            x2 = 0.433395394129247
            x3 = 0.679409568299024
            x4 = 0.865063366688985
            x5 = 0.973906528517172
            w1 = 0.295524224714753
            w2 = 0.269266719309996
            w3 = 0.219086362515982
            w4 = 0.149451349150581
            w5 = 0.066671344308688
            t = [(-x5+1)*0.5, (-x4+1)*0.5, (-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5, (x4+1)*0.5, (x5+1)*0.5]
            w = [w5*0.5, w4*0.5, w3*0.5, w2*0.5, w1*0.5,  w1*0.5, w2*0.5, w3*0.5, w4*0.5, w5*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
        elif n == 11:
            x0 = 0.000000000000000
            x1 = 0.269543155952345
            x2 = 0.519096129110681
            x3 = 0.730152005574049
            x4 = 0.887062599768095
            x5 = 0.978228658146057
            w0 = 0.272925086777901
            w1 = 0.262804544510247
            w2 = 0.233193764591990
            w3 = 0.186290210927734
            w4 = 0.125580369464905
            w5 = 0.055668567116174
            t = [(-x5+1)*0.5, (-x4+1)*0.5, (-x3+1)*0.5, (-x2+1)*0.5, (-x1+1)*0.5, 0.5, (x1+1)*0.5, (x2+1)*0.5, (x3+1)*0.5, (x4+1)*0.5, (x5+1)*0.5]
            w = [w5*0.5, w4*0.5, w3*0.5, w2*0.5, w1*0.5,  w0*0.5, w1*0.5, w2*0.5, w3*0.5, w4*0.5, w5*0.5]
            self.t = np.array(t, dtype = np.float64)
            self.w = np.array(w, dtype = np.float64)
