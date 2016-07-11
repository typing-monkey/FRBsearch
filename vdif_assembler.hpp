#include <iostream>
#include <stdlib.h>

using namespace std;

struct header {
	long unsigned int t0;
	int polarization;
};


struct assembled_chunk {

	long unsigned int t0;
	unsigned char *data;

	assembled_chunk(long int t0);

	~assembled_chunk();

	void set_data(int i, unsigned char x);

};


struct vdif_processor {

	bool upch;
	FILE *output;
	vdif_processor(bool flag1, const char *filename);
	~vdif_processor();
	void process_chunk(assembled_chunk *c);

};



struct vdif_assembler {

	int number_of_processors;
	int start_index, end_index;
	int bufsize;
	int mode;
	string source;
	int port;
	char *filelist_name;
	
	unsigned char *temp_buf;
	unsigned char *data_buf;
	struct header *header_buf;
	queue<assembled_chunk*> chunks;
	
	vdif_processor **processors;
	thread *processor_threads;
	
	vdif_assembler(const char *arg1, const char *arg2);

	~vdif_assembler();

	int register_processor(vdif_processor *p);
	int kill_processor(vdif_processor *p);
	void run();
	void network_capture();
	void read_from_disk();
	void simulate();
	void assemble_chunk();
	assembled_chunk* get_chunk();
	int is_full();
	void vdif_read(unsigned char *data, int size);
	void fill_missing(int n);
};
