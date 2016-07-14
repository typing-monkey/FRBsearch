CC:= g++ -g -mavx -std=c++11 -Wall -Wextra -O3


test: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp upchannelize.cpp
	$(CC) -o test test.cpp -lfftw3f -pthread -lm
