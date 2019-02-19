'''
    Input:
        fname: a string which is the file name to read
    Output:
        If the file contains the right content, returns t0, y0, z0, p0 of the BVP-DAEs
        error: 0 if the file exists
               1 if the file does not exist
'''
def BVPDAEReadData(fname):
    try:
        with open(fname, 'r') as f:
            error = 0
            # read the number of time nodes and the time span variable
            line = f.readline()
            line = f.readline()
            N = int(line)
            line = f.readline()
            line = f.readline()
            list_tspan = line.rstrip().split(' ')
            t0 = np.zeros(N, dtype = np.float64)
            for i in range(N):
                t0[i] = float(list_tspan[i])
            # read the number of y variables and the y variables
            line = f.readline()
            line = f.readline()
            n_y = int(line)
            line = f.readline()
            y0 = np.zeros((N, n_y), dtype = np.float64)
            for i in range(N):
                line = f.readline()
                list_y0 = line.rstrip().split(' ')
                for j in range(n_y):
                    y0[i, j] = float(list_y0[j])
            # read the number of y variables and the y variables
            line = f.readline()
            line = f.readline()
            n_z = int(line)
            line = f.readline()
            z0 = np.zeros((N, n_z), dtype = np.float64)
            for i in range(N):
                line = f.readline()
                list_z0 = line.rstrip().split(' ')
                for j in range(n_z):
                    z0[i, j] = float(list_z0[j])
            # read the number of p variables and the p variables
            line = f.readline()
            line = f.readline()
            n_p = int(line)
            line = f.readline()
            line = f.readline()
            list_p0 = line.rstrip().split(' ')
            p0 = np.zeros(n_p, dtype = np.float64)
            for i in range(n_p):
                p0[i] = float(list_p0[i])
        return error, t0, y0, z0, p0
    except:
        error = 1
        return error

def BVPDAEWriteData(fname, N, n_y, n_z, n_p, tspan, y0, z0, p0):
    with open(fname, 'w') as f:
        # write number of time nodes and the time span in one line
        f.write('{}\n{}\n{}\n'.format('nt:', N, 't:'))
        np.savetxt(f, tspan, delimiter=' ', newline=" ")
        # write number of y variables and the y variables
        f.write('\n{}\n{}\n{}\n'.format('ny:', n_y, 'y:'))
        np.savetxt(f, y0, delimiter=' ')
        # write number of z variables and the z variables
        f.write('{}\n{}\n{}\n'.format('nz:', n_z, 'z:'))
        np.savetxt(f, z0, delimiter=' ')
        # write number of p variables and the p variables in one line
        f.write('{}\n{}\n{}\n'.format('np:', n_p, 'p:'))
        np.savetxt(f, p0, delimiter=' ', newline=" ")