include ./src/rf_pipelines/Makefile.local

CFLAG:= -ggdb -mavx -std=c++11 -O3
LIBCYTHON=vdif_assembler_cython.so

all:test $(LIBCYTHON)

test: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp aro_stream.cpp gaussian.cpp
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lfftw3f -pthread -lrf_pipelines -lm

%_cython.cpp: %_cython.pyx vdif_assembler.hpp
	cython --cplus $<

vdif_assembler_cython.so: vdif_assembler_cython.cpp
	$(CPP) -shared -o $@ $< -lfftw3f -pthread -lrf_pipelines -lm -lhdf5 -lpng 