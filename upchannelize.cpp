<<<<<<< HEAD
#include <iostream>
#include <ipp.h>

void upchannelize(int nt, int nfreq, int size, unsigned char *data, unsigned char *buf) {
	int order = (int)(log((double)size / log(2.0)));
	IppsFFTSpec_C_8sc *pFFTSpec = 0;
	Ipp8u *pFFTSpecBuf, *pFFTInitBuf, *pFFTWorkBuf;
	
	Ipp8sc *pSrc = ippsMalloc_8sc(size);
	Ipp8sc *pDst = ippsMalloc_8sc(size);

	int sizeFFTSpec, sizeFFTInitBuf, sizeFFTWorkBuf;
	ippsFFTGetSize_C_8sc(order, IPP_FFT_DIV_FWD_BY_N, ippAlgHintAccurate, &sizeFFTSpec, &sizeFFTInitBuf, &sizeFFTWorkBuf);

	pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpec);
	pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
	pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

	ippsFFTInit_C_8sc(&pFFTSpec, order, IPP_FFT_DIV_FWD_BY_N, ippAlgHintAccurate, pFFTSpecBuf, pFFTInitBuf);
	if (pFFTInitBuf) ippFree(pFFTInitBuf);

	for (int i = 0; i < (nt/size/2); i++) {

		for (int j = 0; j < nfreq) {
			for (int k = 0; k < size; k++) {
				pSrc[k].re = ((data[i * size * nfreq * 2 + k * 2 * nfreq + j] >> 4) & 0xf) - 8;
				pSrc[k].im = (data[i * size * nfreq * 2 + k * 2 * nfreq + j] & 0xf) - 8;
			}
			ippsFFTFwd_CToC_8sc(pSrc, pDst, pFFTSpec, pFFTWorkBuf);
			for (int k = 0; k < size; k++) {
				*buf = (((pSrc[k].re + 8) & 0xf) << 4) + (pSrc[k].im & 0xf) + 8;
				buf++;
			}
			for (int k = 0; k < size; k++) {
				pSrc[k].re = ((data[i * size * nfreq * 2 + (k * 2 + 1) * nfreq + j] >> 4) & 0xf) -8;
				pSrc[k].im = (data[i * size * nfreq * 2 + (k * 2 + 1) * nfreq + j] & 0xf) - 8;
			}
			ippsFFTFwd_CToC_8sc(pSrc, pDst, pFFTSpec, pFFTWorkBuf);
			for (int k = 0; k < size; k++) {
				*buf = (((pSrc[k].re + 8) & 0xf) << 4) + (pSrc[k].im & 0xf) + 8;
				buf++;
			}
		}
	}

	if (pFFTWorkBuf) ippFree(pFFTWorkBuf);
	if (pFFTSpecBuf) ippFree(pFFTSpecBuf);
	ippFree(pSrc);
	ippFree(pDst);

}
=======
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
>>>>>>> 6911bfda92a9a18f9c5f0d32275c4fefba367ae0
