#include <iostream>
#include <fftw3.h>
#include <cmath>

using namespace std;

void upchannelize(int nt, int nfreq, int size, unsigned char *data, float *buf, fftw_plan &p) {

	///int order = (int)(log((double)size / log(2.0)));
	
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

	for (int i = 0; i < (nt/size/2); i++) {

		for (int j = 0; j < nfreq; j++) {
			for (int k = 0; k < size; k++) {
				in[k][0] = ((data[i * size * nfreq * 2 + k * 2 * nfreq + j] >> 4) & 0xf) - 8;
				in[k][1] = (data[i * size * nfreq * 2 + k * 2 * nfreq + j] & 0xf) - 8;
			}
			fftw_execute_dft(p, in, out);
			for (int k = 0; k < size; k++) {
				*buf = pow(out[k][0],2) + pow(out[k][1],2);
				buf++;
			}
			for (int k = 0; k < size; k++) {
				in[k][0] = ((data[i * size * nfreq * 2 + (k * 2 + 1) * nfreq + j] >> 4) & 0xf) -8;
				in[k][1] = (data[i * size * nfreq * 2 + (k * 2 + 1) * nfreq + j] & 0xf) - 8;
			}
			fftw_execute_dft(p, in, out);
			for (int k = 0; k < size; k++) {
				*buf = pow(out[k][0],2) + pow(out[k][1],2);
				buf++;
			}
		}
	}
	fftw_free(in);fftw_free(out);


}
