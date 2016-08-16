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
	
	void set_time(long unsigned int time);
	void set_data(int i, unsigned char x);

};


struct vdif_processor {
	int p_id;
	bool is_running,upch;
	assembled_chunk *processor_chunk;
	cpu_set_t p_cpuset;	

	int *temp_buf, *intensity_buf;
	float *fft_buf;

	vdif_processor(int id, bool flag);
	~vdif_processor();
	void process_chunk(int *intensity, int index, char &mask);

};

struct vdif_assembler {

	int number_of_processors;
	int i_start_index, i_end_index;
	int bufsize,p_index;

	int mode;
	string source;
	int port;
	char *filelist_name;	
	int chunk_count;
	FILE *output;
	bool upch, write_to_disk;
	//fftwf_complex **in, **out;	
	
	int *intensity_buffer;
	char intensity_buffer_mask;

	int *job_list;
	
	vdif_processor **processors;
	thread *processor_threads;
		
	vdif_assembler(const char *arg1, const char *arg2, bool flag1, bool flag2, int n);

	~vdif_assembler();


	int get_free_processor();
	void run();
	void network_capture();
	void read_from_disk();
	void simulate();
	void assign_chunk();
	void get_intensity_chunk(float *intensity, ssize_t stride);

};
