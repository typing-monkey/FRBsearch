#include <time.h>
#include <sys/socket.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <arpa/inet.h>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <fftw3.h>
#include <rf_pipelines.hpp>

#include "gaussian.cpp"

#include "vdif_assembler.hpp"
#include "square_sum.cpp"
#include "upchannelize.cpp"

using namespace std;
using namespace rf_pipelines;

mutex mtx,mtx2,mtx3;
condition_variable cv,cv_chunk;
fftwf_plan p;
bool upch = false;
bool write_to_disk = false;

char *strptime(const char * __restrict, const char * __restrict, struct tm * __restrict);

namespace constants{

	const int nfreq = 1024;
	const int header_size = 12; //int64 time + int32 thread ID
	const int chunk_size = 65536*8;
	const int max_chunks = 8;
	const int max_processors = 10;
	const int frame_per_second = 390625;
	const unsigned int buffer_size = chunk_size * 3; //at most 3 chunk in buffer
	const int file_packets = 131072;
	const int udp_packets = 32;
	const int nfft = 16;
	const int nsample_integrate = 512;
	const int intensity_chunk_size = (int) chunk_size / nsample_integrate * nfreq;
	const int intensity_buffer_size = intensity_chunk_size * max_chunks;
}


assembled_chunk::assembled_chunk() {
	
	data = new unsigned char[constants::chunk_size * constants::nfreq];
	t0 = 0;

}

assembled_chunk::~assembled_chunk() {
	delete[] data;
}

void assembled_chunk::set_data(int i, unsigned char x) {
	data[i] = x;
}

vdif_processor::vdif_processor(){
	is_running = false;

}

vdif_processor::~vdif_processor(){

}

void vdif_processor::process_chunk(shared_ptr<assembled_chunk> c, int *intensity, fftwf_complex *in, fftwf_complex *out, int index, char &mask) {
	unique_lock<mutex> lk2(mtx2);
	is_running = true;
	lk2.unlock();
	cout << c->t0 << endl;

	unsigned char temp[constants::chunk_size];
	
	for (int i = 0; i < constants::nfreq; i++) {
		
		for (int j = 0; j < constants::chunk_size; j++) {
			temp[j] = c->data[i + constants::nfreq*j];
		}

		int nintegrate = 9;
		if (upch) {
			upchannelize_square_sum(constants::chunk_size, constants::nfft, temp, in, out, p);
			nintegrate -= 4;
		}
		square_sum(temp, constants::chunk_size, nintegrate, intensity + i * (constants::intensity_chunk_size / constants::nfreq));
			
	} 
	
	//gaussian_number(intensity, 1024*1024);

	cout << "Processing chunk done." << endl;

	unique_lock<mutex> lk3(mtx3);
 
	mask = mask | (1 << index);
	
	lk3.unlock();

	lk2.lock();
	is_running = false;
	lk2.unlock();

}


vdif_assembler::vdif_assembler(const char *arg1, const char *arg2, bool flag1, int n){
	number_of_processors = n;
	if (strcmp("network",arg1)==0) {
		temp_buf = new unsigned char[constants::udp_packets * 1056];
		mode = 0;
		port = atoi(arg2);
	} else if (strcmp("disk",arg1)==0) {
		temp_buf = new unsigned char[constants::file_packets * 1056];
		mode = 1;
		filelist_name = new char[strlen(arg2)];
		strcpy(filelist_name,arg2);
	} else if (strcmp("simulate",arg1)==0) {
		temp_buf = new unsigned char[1056];
		mode = 2;
	} else {
		cout << "Unsupported option." << endl;
		exit(1);
	}

	in = new fftwf_complex *[number_of_processors];
	out = new fftwf_complex *[number_of_processors];
	for (int i = 0; i< number_of_processors; i++) {
		in[i] = fftwf_alloc_complex(constants::nfft);
		out[i] = fftwf_alloc_complex(constants::nfft);
	}

	if (flag1) {
		upch = true;		
		p = fftwf_plan_dft_1d(constants::nfft, in[0], out[0], FFTW_FORWARD, FFTW_MEASURE);
	}
	
	number_of_processors = n;
	processors = new vdif_processor *[number_of_processors];
	vdif_processor p[number_of_processors];
	for (int i = 0; i< number_of_processors; i++) {
		processors[i] = &p[i];
	}
	data_buf = new unsigned char[constants::buffer_size * constants::nfreq];
	header_buf = new struct header[constants::buffer_size];
	bufsize = 0;
	chunk_count = 0;
	i_start_index = 0; start_index = 0;
	i_end_index = 0; end_index = 0;
	processor_threads = new thread[number_of_processors];
	intensity_buffer = new int[constants::intensity_buffer_size];
	intensity_buffer_mask = 0;

}

vdif_assembler::~vdif_assembler() {
	
	for (int i = 0; i < number_of_processors; i++) {
		fftwf_free(in[i]);
		fftwf_free(out[i]);
	}

	delete[] processors;
	delete[] data_buf;
	delete[] header_buf;
}


/*
int vdif_assembler::register_processor(vdif_processor *p) {

	if (number_of_processors < constants::max_processors) {

		processors[number_of_processors] = p;
		number_of_processors++;
		return 1;
	} else {
		cout << "The assembler is full, can't register any new processors." << endl;
		return 0;
	}
}

int vdif_assembler::kill_processor(vdif_processor *p) {
	for (int i = 0; i < number_of_processors; i++) {
		if (processors[i] == p) {
			for (int j = i; j < (number_of_processors - 1); j++) {
				processors[j] = processors[j + 1];
			}
			number_of_processors--;
			return 1;
		}
	}
	cout << "Unable to find this processor.";
	return 0;
} */

int vdif_assembler::is_full() {

	return bufsize == constants::buffer_size;
}


void vdif_assembler::run() {

	thread assemble_t(&vdif_assembler::assemble_chunk,this);
	//thread intensity_t(&vdif_assembler::intensity_streamformer,this);
	thread stream_t;

	if (mode==0) {
		stream_t = thread(&vdif_assembler::network_capture,this);
	} else if (mode==1) {
		stream_t = thread(&vdif_assembler::read_from_disk,this);
	} else if (mode==2) {
		stream_t = thread(&vdif_assembler::simulate,this);
	}
	
	stream_t.join();
	//intensity_t.join();
	assemble_t.join();

}

void vdif_assembler::intensity_streamformer() {

	FILE *output;
	output = fopen("test.dat","w");
	if (!output) {
		cout << "Can't write to file." << endl;
		exit(10);
	}

	for (;;) {
		
		unique_lock<mutex> lk3(mtx3);
		if (intensity_buffer_mask & (1 << i_start_index)) {
			fwrite(intensity_buffer + i_start_index * constants::intensity_chunk_size, sizeof(int), constants::intensity_chunk_size, output);
			intensity_buffer_mask = intensity_buffer_mask ^ (1 << i_start_index);
			cout << "Writing chunk " << chunk_count << " to test.dat..." << endl;
			memset(intensity_buffer + i_start_index * constants::intensity_chunk_size, 0 ,constants::intensity_chunk_size * 4);
			i_start_index = (i_start_index + 1) % constants::max_chunks;
			chunk_count++;
		}
		lk3.unlock();
	}
	fclose(output);
}

void vdif_assembler::get_intensity_chunk(int *buf) {

	for (;;) {
		unique_lock<mutex> lk3(mtx3);		
		if (intensity_buffer_mask & (1 << i_start_index)) {

			for (int i = 0; i < constants::intensity_chunk_size; i++ ) {
		
				buf[i] = intensity_buffer[i_start_index * constants::intensity_chunk_size + i];

			}	

			intensity_buffer_mask = intensity_buffer_mask ^ (1 << i_start_index);
			lk3.unlock();
			cout << "Finish assembling chunk " << i_start_index << ", sending to bonsai..." << endl;
			i_start_index = (i_start_index + 1) % constants::max_chunks;
			break;
		}
		lk3.unlock();
	}



}

void vdif_assembler::assemble_chunk() {
	
	for (;;) {

		//cout << " start: " << start_index << " end: " << end_index << endl;
		unique_lock<mutex> lk(mtx);
		
		if (bufsize < constants::chunk_size) {
			cv.wait(lk);
		}
		
		cout << "Chunk found" << endl;
		shared_ptr<assembled_chunk> c = make_shared<assembled_chunk>();
		//cout << header_buf[start_index].t0 << endl;
		c->t0 = header_buf[start_index].t0;	
		bool assigned = false;
		int processor_index = 0;

		for (int i = 0; i < constants::chunk_size * constants::nfreq; i++) {
			c->set_data(i, data_buf[start_index*constants::nfreq+i]);
		}
		
		start_index = (start_index + constants::chunk_size) % constants::buffer_size;
		bufsize -= constants::chunk_size;
		cout << "excess: " << bufsize << endl;

		
		
		while (!assigned) { 
			unique_lock<mutex> lk2(mtx2);
			if (!processors[processor_index]->is_running) {
				
				processor_threads[processor_index] = thread(&vdif_processor::process_chunk,processors[processor_index],c,intensity_buffer+i_end_index*constants::intensity_chunk_size,in[processor_index],out[processor_index],i_end_index,ref(intensity_buffer_mask));
				processor_threads[processor_index].detach();
				assigned = true;
				i_end_index = (i_end_index + 1) % constants::max_chunks;
				
			} else {
				processor_index++;
			}
			lk2.unlock();
			if (processor_index == number_of_processors) {
				cout << "All processors are busy. Waiting for them to finish..." << endl;
				processor_index = 0;
			}

			
		}
		lk.unlock();

		
	}
}


void vdif_assembler::network_capture() {
	
	int size = constants::udp_packets * 1056;
	struct sockaddr_in server_address;
	memset(&server_address, 0, sizeof(server_address));
	
	int sock_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (sock_fd < 0) {
		cout << "socket failed." << std::endl;
	}
	server_address.sin_family = AF_INET;
	server_address.sin_port = htons(port);
	server_address.sin_addr.s_addr = inet_addr("127.0.0.1");
	if (bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address)) < 0) {
		cout << "bind failed." << std::endl;
	}

	for (;;) {
		if (read(sock_fd, temp_buf, size) == size) {
			unique_lock<mutex> lk(mtx);
			if (!is_full()) {
				vdif_read(temp_buf, size);
			}
			else {
				cout << "Buffer is full. Dropping packets." << endl;
			}
			lk.unlock();		
		}
	}

}

void vdif_assembler::read_from_disk() {        
	
	ifstream fl(filelist_name, ifstream::in);
        string filename;
        
        int bytes_read;	
	if (!fl) {
		cout << "Cannot open filelist." << endl;
		exit(1);
	}
        while (getline(fl, filename)){
		bytes_read = 0;
		FILE *fp = fopen(filename.c_str(), "r");
                if (!fp) {
                        cout << "can't open " << filename << endl;
                      	continue; 
                }
		cout << "Reading " << filename << endl;
		
               	while ((bytes_read < constants::file_packets * 1056) && !(feof(fp))) {
			fread(&temp_buf[bytes_read],sizeof(temp_buf[bytes_read]),1,fp);
			bytes_read++;
		}
		
		vdif_read(temp_buf,bytes_read);
                
		fclose(fp);
		
        }
}


void vdif_assembler::simulate() {
	
	int word[8];
	struct tm epoch;
	strptime("2000-01-01 00:00:00","%Y-%m-%dT %H:%M:%S", &epoch);
	
	
	for(int t0 = difftime(time(0),mktime(&epoch));;t0++) {
		//int duration = (int)((rand() % 3 + 3)*1000);
		int duration = 3000;
		//int start_frame = rand() % (constants::frame_per_second - 2000);
		//start_frame = 1000;
		//cout << duration << " " << start_frame << endl;
		for (int frame = 0; frame < constants::frame_per_second; frame++) {
			unsigned char voltage = 136;
			//if ((frame > start_frame) && (frame < start_frame+duration)) {
				//cout << "pulse!" << endl;
			if (((frame % 50000) >= 10000) && ((frame % 50000) < (10000+duration))) {
				voltage = 137;
			}
			for (int pol = 0; pol < 2; pol++) {
				
				word[0] = t0 & 0x3FFFFFFF;
				word[1] = frame & 0xFFFFFF;
				word[3] = (pol << 16) & 0x3FF0000;
				
				for (int i = 0; i < 8; i++){
					for (int j = 0; j < 4; j++) {
						temp_buf[i*4+j] =(int) ((word[i] >> (j*8)) & 0xFF);
					}
				}

					
				for (int i = 32; i < 1056; i++) {
					temp_buf[i] = voltage;
				}
				
				unique_lock<mutex> lk(mtx);
				if (!is_full()) {
					vdif_read(temp_buf, 1056);
				}
				else {
					cout << "Buffer is full. Dropping packets." << endl;
				}
				lk.unlock();
			}	
		}
		this_thread::sleep_for(chrono::seconds(1));
	}

}



void vdif_assembler::vdif_read(unsigned char *data, int size) {


	int word[8];
	int count = 0;
	long int t0 = 0;
	int pol = 0;
	int nmissing = 0;
	long int current,expect = 0;
	bool invalid;
	
	while ((count < size) && (!is_full())) {
		
		for (int i = 0; i < 8; i++) {
			word[i] = (data[count + 3] << 24) + (data[count + 2] << 16) + (data[count + 1] << 8) + data[count];
			count += 4;
		}
		
		invalid = (word[0] >> 31);
		
 		if (invalid) {
			fill_missing(1);
		}
		
		t0 = (long int) (word[0] & 0x3FFFFFFF) * (long int) constants::frame_per_second + (long int) (word[1] & 0xFFFFFF);
		pol = (word[3] >> 16) & 0x3FF;
		
		current = t0 * 2 + pol;

		if (expect) {
			nmissing = (int) (current - expect);
		}

		if (nmissing) {
			cout << "current: " << current << " expect: " << expect << endl; 
			cout << "start: " << start_index << " end: " << end_index << endl;
			fill_missing(nmissing);
		}

		expect = current + 1;
		
		header_buf[end_index].t0 = t0;
		header_buf[end_index].polarization = pol;
		
		for (int i = 0; i < constants::nfreq; i++) {
			data_buf[end_index * constants::nfreq + i] = data[count];
			count++;
		}
		end_index = (end_index + 1) % constants::buffer_size;
		bufsize++;
	
		if (bufsize >= constants::chunk_size) {
			cv.notify_one();
		}

		//cout << "start: " << start_index << " end: " << end_index << " size: " << bufsize << endl;
	}
}

void vdif_assembler::fill_missing(int n) {
	long int prev_t0 = header_buf[(end_index-1) % constants::buffer_size].t0;
	int prev_pol = header_buf[(end_index-1) % constants::buffer_size].polarization;
	//cout << "Missing " << n << " packets." << endl;
	for (int i = 0; i < n; i++) {
		prev_t0++;
		prev_pol = (prev_pol + 1) % 2;
		header_buf[end_index].t0 = prev_t0;
		header_buf[end_index].polarization = prev_pol;
		for (int j = 0; j < constants::nfreq; j++) {
			data_buf[end_index * constants::nfreq + j] = 0;
		}
		end_index = (end_index + 1) % constants::buffer_size;
		bufsize++;
	}
	
}


	
