#include <fftw3.h>
#include <iostream>
#include <math.h>

using namespace std;

void upchannelize_square_sum(int nsample, int N, unsigned char *data,fftwf_complex *in, fftwf_complex *out, fftwf_plan &p) {
	
	if (nsample % N) {
		cout << "Number of samples to FFT must evenly divide FFT size." << endl;
		exit(10);
	}
	
	//p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nsample / N; i++) {
		for (int j = 0; j < N; j++) {
			in[j][0] = float((data[i * 16 + j] >> 4) - 8);
			in[j][1] = float((data[i * 16 + j] & 0xF) - 8);
			if (isnan(in[j][0])) {
				cout << in[j][0] << endl;
			}
			if (isnan(in[j][1])) {
				cout << in[j][1] << endl;
			}
		}
		fftwf_execute_dft(p, in, out);
		for (int j = 0; j < N; j++) {
			data[nsample / N*j + i] = (unsigned char)((int(round(out[j][0] / N)) << 4) + (int(round(out[j][1] / N))));
		}
	}

	//fftwf_destroy_plan(p);
}
