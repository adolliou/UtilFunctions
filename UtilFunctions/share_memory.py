import numpy as np
from multiprocessing.shared_memory import SharedMemory


def gen_shmm(create=False, name=None, ndarray=None, size=0, shape=None, dtype=float):
    assert (type(ndarray) != type(None) or size != 0) or type(name) != type(None)
    assert type(ndarray) != type(None) or type(shape) != type(None)
    size = size if type(ndarray) == type(None) else ndarray.nbytes
    shmm = SharedMemory(create=create, size=size, name=name)
    shmm_data = np.ndarray(shape=shape if type(ndarray) == type(None) else ndarray.shape,
                           buffer=shmm.buf, dtype=dtype)

    if create and type(ndarray) != type(None):
        shmm_data[:] = ndarray[:]
    elif create:
        shmm_data[:] = np.nan

    return shmm, shmm_data

# example using shared variable:
# self._shmm_par,self.data_par = gen_shmm(create=True,ndarray=self.data_par)
# self._par = {"name":self._shmm_par.name ,"dtype":self.data_par.dtype ,"shape":self.data_par.shape }
# shmm_war,data_war = gen_shmm(create=False,**war)
