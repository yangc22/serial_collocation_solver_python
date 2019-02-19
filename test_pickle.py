import pickle
import numpy as np
import pprint

y = np.zeros((100, 5))
size_y = 5

# with open("test1.data", 'wb') as data:
with open('test.data', 'wb') as data:
    size_y_str = 'ny = \n'
    y_str = 'y = \n'
    pickle.dump(size_y_str, data, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(size_y, data, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(y_str, data, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(y, data, protocol=pickle.HIGHEST_PROTOCOL)
    # data.write(size_y_str)
    # data.write('%d' %(size_y))
    # data.write(y_str)
    # np.savetxt(data, y, fmt='%-7.2f')
    # np.savetxt(data, y, fmt='%.5f')
    np.savetxt(data, y)
data.close()
# np.savetxt("test.txt", y)
'''
read_file = open('test.data', 'rb')
while True:
        try:
            pprint.pprint(pickle.load(read_file))
        except EOFError:
            break
'''
'''
with open('test.data', 'rb') as file:
    try:
        while True:
            yield pickle.load(file)
    except EOFError:
        pass
'''