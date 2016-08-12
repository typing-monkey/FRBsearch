include ./src/rf_pipelines/Makefile.local

CFLAG:= -ggdb -mavx -std=c++11 -O3

all:test

test: test.cpp vdif_assembler.hpp vdif_assembler.cpp aro_stream.cpp upchannelize.cpp gaussian.cpp preprocessing.cpp
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lfftw3 -pthread -lrf_pipelines -lm


