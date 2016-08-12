include Makefile.local

#CC:= g++ -g -Wall -Wextra -Wconversion

INCFILES= vdif_assembler.hpp
#CFLAG:= -ggdb -mavx -std=c++11 -O3
LIBCYTHON=vdif_assembler_cython.so
PYMODULES=vdif_assembler.py
#OFILES=vdif_assembler.o
LIBFILES=libvdif_assembler.so
BINFILES=test

#all: $(BINFILES) $(LIBFILES) $(LIBCYTHON) $(TESTBINFILES)
all: libvdif_assembler.so vdif_assembler_cython.so

#cython: $(LIBCYTHON)

# %.o: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp aro_stream.cpp gaussian.cpp
	# $(CPP) $(CPP_LFLAGS) -c -o $@ $< -I$(PYINCDIR) -lfftw3f -pthread -lrf_pipelines -lpng -lm

# %.o: %.cpp $(INCFILES)
# 	$(CPP) -c -o $@ $< -lrf_pipelines -lm

# vdif_assembler.o: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp aro_stream.cpp gaussian.cpp
# 	$(CPP) -c -o $@ $< -lrf_pipelines -lm

# test: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp aro_stream.cpp gaussian.cpp preprocessing.cpp
# 	$(CPP) $(CPP_LFLAGS) -c -o $@ $< -I$(PYINCDIR) -lfftw3f -pthread -lrf_pipelines -lpng -lm

vdif_assembler.o: vdif_assembler.cpp vdif_assembler.hpp
	$(CPP) -c -o $@ $<
	#$(CPP) $(CPP_LFLAGS) -c -o $@ $< -I$(PYINCDIR) -lfftw3f -pthread -lrf_pipelines -lpng -lm

vdif_assembler_cython.cpp: vdif_assembler_cython.pyx vdif_assembler_pxd.pxd vdif_assembler_cython.hpp $(INCFILES)
	cython --cplus $<

libvdif_assembler.so: vdif_assembler.o
	$(CPP) -o $@ -shared $^

vdif_assembler_cython.so: vdif_assembler_cython.cpp libvdif_assembler.so vdif_assembler.hpp
	$(CPP) -shared -o $@ $<  -lvdif_assembler -lfftw3f -pthread -lrf_pipelines -lhdf5 -lpng -lm 

install: $(INCFILES) $(LIBFILES) $(LIBCYTHON)
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(LIBFILES) $(LIBDIR)/
	cp -f $(LIBCYTHON) $(PYMODULES) $(PYDIR)/
	#cp -f vdif_assembler_cython.cpp $(PYMODULES) $(PYDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(LIBFILES); do rm -f $(LIBDIR)/$$f; done
	for f in $(LIBCYTHON) $(PYMODULES); do rm -f $(PYDIR)/$$f; done

clean:
	rm -f *~ *.o *_cython.cpp *.so *.pyc $(BINFILES) $(TESTBINFILES)

# include ./src/rf_pipelines/Makefile.local

# CFLAG:= -ggdb -mavx -std=c++11 -O3
# IPP_DIR = /home/surp2016lyu/ipp/include

# all:test

# test: test.cpp vdif_assembler.hpp vdif_assembler.cpp aro_stream.cpp gaussian.cpp preprocessing.cpp
# 	$(CPP) $(CPP_LFLAGS) -o $@ $< -I$(IPP_DIR) -lfftw3 -pthread -lrf_pipelines -lm

