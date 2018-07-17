from collocation_coefficients import *

for i in [3,4,5,6,7,8,9,10]:
    rk = lobatto(i)
    print(rk.A.shape)