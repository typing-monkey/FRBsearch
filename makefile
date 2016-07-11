CC:= g++ -mavx -std=c++11 -Wall -Wextra -O3


test: test.cpp vdif_assembler.hpp vdif_assembler.cpp square_sum.cpp
	$(CC) -o test test.cpp -I "/usr/local/lib" -L "/usr/local/include"