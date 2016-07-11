#include <fftw3.h>
#include <iostream>
#include <math.h>

using namespace std;

void upchannelize_square_sum(int nsample, int N, unsigned char *data) {
	
	if (nsample % N) {
		cout << "Number of samples to FFT must evenly divide FFT size." << endl;
		exit(10);
	}
	
	fftwf_complex *in, *out;
	fftwf_plan p;
	in = fftwf_alloc_complex(N);
	out = fftwf_alloc_complex(N);

	p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nsample / N; i++) {
		for (int j = 0; j < N; j++) {
			in[j][0] = float((data[i * 16 + j] >> 4) - 8);
			in[j][1] = float((data[i * 16 + j] & 0xF) - 8);
		}
		fftwf_execute(p);
		for (int j = 0; j < N; j++) {
			data[nsample / N*j + i] = (unsigned char)((int(round(out[j][0] / N)) << 4) + (int(round(out[j][1] / N))));
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(in); fftwf_free(out);
}