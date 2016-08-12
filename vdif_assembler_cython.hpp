#include <vdif_assembler.hpp>

struct cython_assembler{
	vdif_assembler *a;
	cython_assembler(const char* arg1, const char* arg2, bool flag, bool todisk, int n){
		a = new vdif_assembler(arg1,arg2,flag,todisk,n);
	}

	void run(){
		a->run();
	}

	void get_intensity_chunk(float* buf,int ntime){
		a->get_intensity_chunk(buf,ntime);
	}
};