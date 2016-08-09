from libcpp cimport bool

cdef extern from "vdif_assembler_cython.hpp":
	cpdef cppclass cython_assembler:
		cython_assembler(const char *arg1, const char *arg2, bool flag1, bool flag2, int n) except +
		void run() except +
		void get_intensity_chunk(int* buf) except +