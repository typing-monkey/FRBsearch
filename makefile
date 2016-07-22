include ./src/rf_pipelines/Makefile.local

CFLAG:= -g -mavx -std=c++11 -O3

all:test

test: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp aro_stream.cpp gaussian.cpp
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lfftw3f -pthread -lrf_pipelines -lm


