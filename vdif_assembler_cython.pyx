import numpy as np
import time
cimport numpy as np

cimport vdif_assembler_pxd
from libcpp cimport bool

dt_native = 2.56*(10**(-6))
n_integrate = 384
n_chunk = 1024
dt_chunk = n_integrate*n_chunk*dt_native

cdef class assembler:
	first_t0 = None
	cdef vdif_assembler_pxd.cython_assembler *p

	def __cinit__(self, char* arg1, char* arg2, bool flag, int n):
		self.p = new vdif_assembler_pxd.cython_assembler(arg1,arg2,flag,<bool>False,n)

	def __dealloc__(self):
		if self.p != NULL:
			del self.p
			self.p = NULL

	def get_intensity_chunk(self):
		#assert empty_ar.dtype == np.int32_t, "get_intensity_chunk got incorrect dtype array"
		cdef np.ndarray[np.int32_t, ndim=2] data_ar = np.empty((1024,1024),dtype=np.int32)
		self.p.get_intensity_chunk(<int*> data_ar.data)
		if self.first_t0 == None:
			self.first_t0 = time.time()

		self.first_t0 += dt_chunk
		return (self.first_t0 - dt_chunk), data_ar

	def run(self):
		self.p.run()
