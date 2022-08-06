from netCDF4 import Dataset
import numpy as np
np.set_printoptions(suppress=False)

file_path_original = '/Users/murali/phd/one_phonon_raman/si/bse/si_data/nscf_wo/raman/ndb.BS_elph'

file_path_test = 'nc.temp'

def get_data(filename):
    db = Dataset(filename)
    db.set_always_mask(False)
    db = db['exc_elph'][...]
    return db[...,0] + 1j*db[...,1]


original = get_data(file_path_original)
test = get_data(file_path_test)

print(np.allclose(test,original))
print(test.ravel().shape[0])
print("{:e}".format(test[0,2,999,999]))
print("{:e}".format(test[0,1,348,749]))

reshaped = test.reshape(6,500,2,125,4)
print('Reshape: {:e}'.format(reshaped[3,127,1,73,2]))

transposed = reshaped.transpose(4,2,0,1,3)
print('Transpose: {:e}'.format(transposed[2,1,3,273,39]))


sliced_arr = transposed[2:4,0,3:5,179:500:3,53:117:5]
print('Sliced: {:e}'.format(sliced_arr[1,0,37,7]))

stripped_arr = sliced_arr[1,0,:,:]
print('Stripped: {:e}'.format(stripped_arr[13,5]))

matmul = test[0,1,:,:].T@np.conj(test[0,2,:,:].T)
einsum= np.einsum("ijkk,ijlk->ijl",test,test)

#print('Sliced: {:e}'.format())

print('matmul: {:e}'.format(matmul[239,739]))
print('einsum: {:e}'.format(einsum[0,2,473]))