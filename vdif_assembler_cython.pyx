import numpy as np
cimport numpy as np

cimport vdif_assembler_pxd
from libcpp cimport bool

cdef class assembler:
	cdef vdif_assembler_pxd.cython_assembler *p

	def __cinit__(self, char* arg1, char* arg2, bool flag, int n):
		self.p = new vdif_assembler_pxd.cython_assembler(arg1,arg2,flag,n)

	def __dealloc__(self):
		if self.p != NULL:
			del self.p
			self.p = NULL

	def get_intensity_chunk(self, np.ndarray[np.int32_t, ndim=2] empty_ar):
		#assert empty_ar.dtype == np.int32_t, "get_intensity_chunk got incorrect dtype array"
		self.p.get_intensity_chunk(<int*> empty_ar.data)

	def run(self):
		self.p.run()
