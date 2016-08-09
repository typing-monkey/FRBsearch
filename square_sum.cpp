#include <stdio.h>
#include <iostream>

#include "ch_vdif_assembler_kernels.hpp"

using namespace std;

void square_sum(unsigned char *raw_data, int nt, int nsample_integrate, int *intensity) {

	if (nt % (1 << nsample_integrate)) {
		cout << "Number of sampels to integrate must evenly divide number of samples." << endl;
		exit(10);
	}

	if (nsample_integrate < 4) {
		cout << "Number of samples to sum must evenly divide by 16." << endl;
	}
	int count = 0;
	int sum = 0;
	for (int i = 0; i < nt; i+=16) {
		ch_vdif_assembler::_sum16_auto_correlations(sum, count, &(raw_data[i]));
		intensity[i>>nsample_integrate] += sum;
	}
}
