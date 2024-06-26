import numpy as np
from multiprocessing.shared_memory import SharedMemory


def gen_shmm(create=False, name=None, ndarray=None, size=0, shape=None, dtype=None):
    assert (type(ndarray) != type(None) or size != 0) or type(name) != type(None)
    assert type(ndarray) != type(None) or type(shape) != type(None)
    dtype = ndarray.dtype if dtype is None else dtype = dtype
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

# Processes = []
# for i in range(blabla):
#     Processes.append(Process(target=self.task_fit_pixel,kwargs=keywords))
#     Processes[i].start()
#     if self.verbose>=1: print(f"Starting process job: {i+1:02d} on raster fits\nJob list contains: {len(self.Job_index_list[i])} pixel")
# for i in range(self.Jobs): #join all processes
#     Processes[i].join()

# (inside process function: can use lock to ensure no other process is running while writing on shared memory
# shmm_par, data_par = gen_shmm(create=False, **par)

# lock.acquire()
# data_par[  :,0,i_y,i_x] = best_par #the result UUUUUUgh finally it's here every pixel will be here
# data_cov[:,:,0,i_y,i_x] = best_cov #the result UUUUUUgh finally it's here every pixel will be here
# data_con[    0,i_y,i_x] = best_con #the result UUUUUUgh finally it's here every pixel will be here
# lock.release()