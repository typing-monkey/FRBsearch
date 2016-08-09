#include "rf_pipelines_internals.hpp"
#include "vdif_assembler.cpp"

using namespace std;
condition_variable cv_stream;
mutex mtx_stream;

namespace rf_pipelines {


class aro_stream: public wi_stream {
	
	int *buf;
	vdif_assembler *assembler;

public:
	aro_stream(int nfreq_, double freq_lo_MHz_, double freq_hi_MHz_, double dt_sample_, const char *arg1, const char *arg2, bool flag1, bool flag2, int n) {
		
		this->nfreq = nfreq_;
		this->freq_lo_MHz = freq_lo_MHz_;
		this->freq_hi_MHz = freq_hi_MHz_;
		this->dt_sample = dt_sample_;	
		this->nt_maxwrite = 1024;
		assembler = new vdif_assembler(arg1, arg2, flag1, flag2, n);

		buf = new int[nfreq * nt_maxwrite];
	}
	
	~aro_stream() {
		delete[] buf;
	}
	
	virtual void stream_body(wi_run_state &run_state) override {
		
		run_state.start_substream(0.0);
		thread main_t(&vdif_assembler::run, assembler);
					
		for (;;) {	
			float *intensity;
			float *weights;
			ssize_t stride;
			bool zero_flag = false;
			
			run_state.setup_write(nt_maxwrite, intensity, weights, stride, zero_flag);
			
			assembler->get_intensity_chunk(buf);
		
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < nt_maxwrite; j++) {		
					
					intensity[i*stride + j] = buf[i*nfreq+j];
					weights[i*stride + j] = 1.0;
					
				}
			}
			//cout << "Bonsai received a chunk." << endl;

			run_state.finalize_write(nt_maxwrite);
					
		}
		run_state.end_substream();
		
		main_t.join();
	}

};


}
