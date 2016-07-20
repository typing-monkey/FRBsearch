#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace rf_pipelines;

struct header {
	long unsigned int t0;
	int polarization;
};


struct assembled_chunk {
	long unsigned int t0;
	unsigned char *data;

	assembled_chunk();

	~assembled_chunk();

	void set_data(int i, unsigned char x);

};


struct vdif_processor {

	bool is_running;
	vdif_processor();
	~vdif_processor();
	void process_chunk(shared_ptr<assembled_chunk> c, int *intensity, fftwf_complex *in, fftwf_complex *out, int index, char &mask);

};



struct vdif_assembler {

	int number_of_processors;
	int start_index, i_start_index, i_end_index, end_index;
	int bufsize,i_bufsize;
	int mode;
	string source;
	int port;
	char *filelist_name;	
	int chunk_count;

	fftwf_complex **in, **out;	

	unsigned char *temp_buf;
	unsigned char *data_buf;
	struct header *header_buf;
	
	int *intensity_buffer;
	char intensity_buffer_mask;
	
	vdif_processor **processors;
	thread *processor_threads;
		
	vdif_assembler(const char *arg1, const char *arg2, bool flag1, int n);

	~vdif_assembler();

	//int register_processor(vdif_processor *p);
	//int kill_processor(vdif_processor *p);
	void run();
	void intensity_streamformer();
	void network_capture();
	void read_from_disk();
	void simulate();
	void assemble_chunk();
	void get_intensity_chunk(int *buf);
	int is_full();
	void vdif_read(unsigned char *data, int size);
	void fill_missing(int n);
};
