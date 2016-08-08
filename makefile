include ./src/rf_pipelines/Makefile.local

CFLAG:= -ggdb -mavx -std=c++11 -O3
IPP_DIR = /home/surp2016lyu/ipp/include

all:test

test: test.cpp vdif_assembler.hpp vdif_assembler.cpp aro_stream.cpp gaussian.cpp preprocessing.cpp
	$(CPP) $(CPP_LFLAGS) -o $@ $< -I$(IPP_DIR) -lfftw3 -pthread -lrf_pipelines -lm


