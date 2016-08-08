#include <fftw3.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <png++/png.hpp>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;

const float pi = (atan(1)*4.0);

float blackman_fun(int index, int width){
	float coef = static_cast<float>(index)/static_cast<float>(width);
	return 0.42 - 0.5*cos(2.0*pi*coef) + 0.08*cos(4*pi*coef);
}

void neg_normalize(float* ar, ssize_t len){
	float sumf = 0;
	for(ssize_t i = 0; i < len; i++){
		sumf += ar[i];
	}
	for(ssize_t i = 0; i < len; i++){
		ar[i] = -ar[i]/sumf;
	}
}

void mult_inplace(float* a, float* b, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i] = a[i]*b[i];
	}
}

void sub_inplace(float* a, float* b, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i] -= b[i];
	}
}

void sqrt_inplace(float* a, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i] = sqrt(a[i]);
	}
}

void add_inplace(float* a, float* b, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i] += b[i];
	}
}

void div_inplace(float* a, float* b, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i] /= b[i];
	}
}

float sum(float* ar, ssize_t len){
	float sumf = 0;
	for(ssize_t i = 0; i < len; i++){
		sumf += ar[i];
	}
	return sumf;
}

float mean(float* ar, ssize_t len){
	return sum(ar,len)/((float) len);
}

float* broadcast(float val, ssize_t len){
	float* ret = new float[len];
	for(ssize_t i = 0; i < len; i++){
		ret[i] = val;
	}
	return ret;
}

float* mult(float* a, float* b, ssize_t len){
	float* ret = new float[len];
	for(ssize_t i = 0; i < len; i++){
		ret[i] = a[i]*b[i];
	}
	return ret;
}

float* add(float* a, float* b, ssize_t len){
	float* ret = new float[len];
	for(ssize_t i = 0; i < len; i++){
		ret[i] = a[i] + b[i];
	}
	return ret;
}

float sample_var(float* ar, ssize_t len){
	float flen = (float) len;
	float* interim_sum = mult(ar,ar,len);
	float ret = ((flen - 1)/flen)*sum(interim_sum,len);
	delete[] interim_sum;
	return ret;
}

float skew(float*ar, ssize_t len){
	float flen = (float) len;
	float* inter_sum = mult(ar,ar,len);
	mult_inplace(inter_sum,ar,len);
	float ret = ((flen - 1)/flen)*sum(inter_sum,len);
	delete[] inter_sum;
	return ret;
}


float median(float* in, ssize_t len){
	vector<float> t_vec(in, in + len);
	sort(t_vec.begin(),t_vec.end());
	if(len % 2 == 0){
		ssize_t mid = len/2;
		return 0.5*(t_vec[mid-1] + t_vec[mid]);
	}
	else{
		return t_vec[(len - 1)/2];
	}
}

void mult_inplace_complex(fftw_complex* a, fftw_complex* b, ssize_t len){
	for(ssize_t i = 0; i < len; i++){
		a[i][0] = a[i][0]*b[i][0] - a[i][1]*b[i][1];
		a[i][1] = a[i][0]*b[i][1] + a[i][1]*b[i][0];
	}
}

void do_stdskew_mask(ssize_t nf, ssize_t nt, float sigma_cut, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	float* vars = new float[nf];
	float* skews = new float[nf];
	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		float* a = mult(&intensity[i*stride],&weights[i],nt);
		float* b = mult(&intensity[i*stride],&weights[i],nt);
		vars[i] = sample_var(a,nt);
		skews[i] = skew(b,nt);
		delete[] a;
		delete[] b;
	}


	// uint8_t* chan_mask = new uint8_t[nf];
	// uint8_t* chan_mask_skew = new  uint8_t[nf];
	float meanvar =  mean(vars,nf);
	float stdvar = sqrt(sample_var(vars,nf));
	float meanskew =  mean(skews,nf);
	float stdskew = sqrt(sample_var(skews,nf));
	int n_rep = 0;
	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		if(vars[i] > sigma_cut*stdvar + meanvar){
			weights[i] = 0;
			n_rep++;
		}

		if(skews[i] > sigma_cut*stdskew + meanskew){
			weights[i] = 0;
			n_rep++;
		}
	}
	delete[] vars;
	delete[] skews;
	cout << "replaced " << n_rep << " channels!" << endl;
}

struct remove_noisy_freq : public wi_transform{
	float sigma_cut;

	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	//ssize_t nt_chunk;
	remove_noisy_freq(const float sigma_cut);
	~remove_noisy_freq(){}
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

remove_noisy_freq::remove_noisy_freq(const float my_sigma_cut){
	sigma_cut = my_sigma_cut;
}

void remove_noisy_freq::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void remove_noisy_freq::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	for(ssize_t i = 0; i < 2; i++){
		do_stdskew_mask(this->nfreq,this->nt_chunk,sigma_cut,intensity,weights,stride,pp_intensity,pp_weight,pp_stride);
		//cout << "ok" << endl;
	}
}

struct remove_continuum : public wi_transform{
	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	//ssize_t nt_chunk;
	remove_continuum(){}
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

void remove_continuum::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void remove_continuum::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	//mean-sub and 
	//sum continuum
	ssize_t nf = this->nfreq;
	ssize_t nt = this->nt_chunk;
	float* freq_sum = new float[nt] ();
	for(ssize_t i = 0; i < nf; i++){
		float ch_mean = mean(&intensity[i*stride],nt);
		float* bc = broadcast(ch_mean,nt);
		sub_inplace(&intensity[i*stride],bc,nt);
		add_inplace(freq_sum,&intensity[i*stride],nt);
		delete[] bc;
	}

	float* mul = mult(freq_sum,freq_sum,nt);
	float* bc = broadcast(sqrt(sum(mul,nt)),nt);
	div_inplace(freq_sum,bc,nt);
	delete[] mul;
	delete[] bc;

	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		float* mul = mult(&intensity[i*stride],freq_sum,nt);
		float* bc = broadcast(sum(mul,nt),nt);
		float* mul2 = mult(bc,freq_sum,nt);
		sub_inplace(&intensity[i*stride],mul2,nt);
		delete[] mul;
		delete[] bc;
		delete[] mul2;
	}
	delete[] freq_sum;
}

struct remove_outliers : public wi_transform{
	float sigma_cut;

	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	//ssize_t nt_chunk;
	remove_outliers(const float my_sigma_cut);
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

remove_outliers::remove_outliers(const float my_sigma_cut){
	sigma_cut = my_sigma_cut;
}

void remove_outliers::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void remove_outliers::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	ssize_t nf = this->nfreq;
	ssize_t nt = this->nt_chunk;
	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		float ch_mean = mean(&intensity[i*stride],nt);
		float ch_std = sqrt(sample_var(&intensity[i*stride],nt));
		for(ssize_t j = 0; j < nt; j++){
			if(abs(intensity[i*stride + j] - ch_mean) > sigma_cut*ch_std){
				intensity[i*stride + j] = ch_mean;
			}
		}
	}
	//fin
}

struct sys_temperature_bandpass : public wi_transform{
	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	//ssize_t nt_chunk;
	sys_temperature_bandpass();
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

sys_temperature_bandpass::sys_temperature_bandpass(){}

void sys_temperature_bandpass::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void sys_temperature_bandpass::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	ssize_t nf = this->nfreq;
	ssize_t nt = this->nt_chunk;
	float* meds = new float[nf];
	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		meds[i] = median(&intensity[i*stride],nt);
	}
	
	float medmed = median(meds,nf);
	//float* zero_ar = broadcast(0, nt);
	#pragma omp parallel for
	for(ssize_t i = 0; i < nf; i++){
		if(meds[i] < medmed*0.001){
			memset(&intensity[i*stride],0,nt);
		}
	}
	delete[] meds;
	//delete zero_ar;
	//fin
}

struct remove_bad_times : public wi_transform{
	float sigma_cut;

	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	remove_bad_times(const float my_sigma_cut);
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

remove_bad_times::remove_bad_times(const float my_sigma_cut){
	sigma_cut = my_sigma_cut;
}

void remove_bad_times::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void remove_bad_times::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	ssize_t nt = this->nt_chunk;
	ssize_t nf = this->nfreq;
	float* stds = new float[nt];
	float* means = new float[nf];
	for(ssize_t i = 0; i < nt; i++){
		float* this_t = new float[nt];
		for(ssize_t j = 0; j < nf; j++){
			this_t[j] = intensity[j*stride + i];
		}
		stds[i] = sqrt(sample_var(this_t,nf));
		delete[] this_t;
	}

	for(ssize_t i = 0; i < nf; i++){
		means[i] = mean(&intensity[i*stride],nt);
	}

	float stdstd = sqrt(sample_var(stds,nt));
	float meanstd = mean(stds,nt);
	int reps = 0;
	for(ssize_t i = 0; i < nt; i++){
		if(stds[i] - meanstd > sigma_cut*stdstd){
			reps++;
			for(ssize_t j = 0; j < nf; j++){
				intensity[j*stride + i] = means[j];
			}
		}
	}
	cout << "replaced " << reps << " time channels!" << endl;

	delete[] stds;
	delete[] means;
}
/*
struct png_writer : public wi_transform{
	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	int chunk = 0;

	char* name_base;
	png_writer(char* name);
	~png_writer(){}
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};
png_writer::png_writer(char* my_name_base){
	name_base = my_name_base;
}

void png_writer::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
}

void png_writer::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	png::image< png::rgb_pixel > image(128, 128);


	//char name[]
	strcat(name_base,to_string(chunk));
	strcat(name_base,".png");
	image.write("rgb.png");
	chunk++;
}
*/


struct highpass_filter : public wi_transform{
	ssize_t n;
	fftw_complex *window_comp, *window_fft, *nt_input, *nt_fft;
	fftw_plan forward_p;
	fftw_plan reverse_p;

	int window_width;
	double freq_lo_MHz;
	double freq_hi_MHz;
	double dt_sample;
	ssize_t nt_maxwrite;
	//ssize_t nt_chunk;
	highpass_filter(const int width);
	~highpass_filter();
	virtual void set_stream(const wi_stream &stream);
	virtual void start_substream(int isubstream, double t0){}
	virtual void end_substream(){}
	virtual void process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride);
};

highpass_filter::highpass_filter(const int width){
	window_width = static_cast<int>(static_cast<float>(width)/0.4054785f);
	if(window_width % 2 == 1){
		window_width++;
	}
}

highpass_filter::~highpass_filter(){
	delete[] window_comp;
	delete[] window_fft;
	delete[] nt_input;
	delete[] nt_fft;
	//delete[] forward_p;
	//delete[] reverse_p;
}

void highpass_filter::set_stream(const wi_stream &stream){
	this->nfreq = stream.nfreq;
	this->freq_lo_MHz = stream.freq_lo_MHz;
	this->freq_hi_MHz = stream.freq_hi_MHz;
	this->dt_sample = stream.dt_sample;
	this->nt_maxwrite = stream.nt_maxwrite;
	this->nt_chunk = stream.nt_maxwrite;
	window_width = min(window_width, (int) this->nt_maxwrite);
	this->nt_prepad = window_width;
	
	n = this->nt_maxwrite + this->nt_prepad;

	window_comp = fftw_alloc_complex(n);
	window_fft = fftw_alloc_complex(n);
	nt_input = fftw_alloc_complex(n);
	nt_fft = fftw_alloc_complex(n);
	fftw_plan window_p = fftw_plan_dft_1d(n, window_comp, window_fft, FFTW_FORWARD, FFTW_ESTIMATE);
	forward_p = fftw_plan_dft_1d(n, nt_input, nt_fft, FFTW_FORWARD, FFTW_MEASURE);
	reverse_p = fftw_plan_dft_1d(n, nt_fft, nt_input, FFTW_BACKWARD, FFTW_MEASURE);

	float* window = new float[n] ();
	for(ssize_t i = 0; i < window_width; i++){
		window[i] = blackman_fun(i,window_width);
	}

	neg_normalize(window,window_width);
	window[window_width/2] += 1;

	//copy window into buffer
	for(ssize_t i = 0; i < n; i++){
		window_comp[i][0] = window[i];
		window_comp[i][1] = 0;
	}

	fftw_execute(window_p);
	//delete[] window_p;
	delete[] window;

	cout << "highpass filter configured with nt_maxwrite " << stream.nt_maxwrite << endl;
	cout << "sample delta t: " << stream.dt_sample << endl;
	// this->nt_prepad = stream.nt_prepad;
	// this->nt_postpad = stream.nt_postpad;
}

void highpass_filter::process_chunk(double t0, double t1, float* intensity, float* weights, ssize_t stride, float* pp_intensity, float* pp_weight, ssize_t pp_stride){
	//cout << "someone called me!" << endl;
	//cout << "found " << (t1 - t0) << " seconds of data" << endl;
	//cout << pp_stride << " " << stride << endl;
	ssize_t nt_p = this->nt_prepad;
	for(ssize_t i = 0; i < this->nfreq; i++){
		for(ssize_t j = 0; j < nt_prepad; j++){
			nt_input[j][0] = pp_intensity[i*pp_stride + j];
			nt_input[j][1] = 0;

			nt_fft[j][0] = 0;
			nt_fft[j][1] = 0;
		}
		for(ssize_t j = 0; j < this->nt_chunk; j++){
			nt_input[j + nt_p][0] = intensity[i*stride + j];
			nt_input[j + nt_p][1] = 0;

			nt_fft[j + nt_p][0] = 0;
			nt_fft[j + nt_p][1] = 0;
		}
		fftw_execute(forward_p);
		mult_inplace_complex(nt_fft, window_fft, n);
		fftw_execute(reverse_p);
		for(ssize_t j = 0; j < stride; j++){
			intensity[i*stride + j] = nt_input[j + nt_p][0];
		}
	}

	//cout << "highpass filter finished" << endl;
}
