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
#include <rf_pipelines.hpp>
#include <sys/time.h>
#include <chrono>
//#include "gaussian.cpp"
#include <fftw3.h>

#include <sched.h>
#include <pthread.h>

#include "vdif_assembler.hpp"
#include "square_sum_andre.cpp"
//#include "upchannelize.cpp"
#include <assert.h>


using namespace std;
using namespace rf_pipelines;

//global variable declaration
auto start = chrono::steady_clock::now();
mutex mtx,mtx2,mtx3;
condition_variable *cv;
float perf, tot_perf, avg_perf; //performance measurer
fftw_plan p;

//system time
char *strptime(const char * __restrict, const char * __restrict, struct tm * __restrict);

//constants are defined here
namespace constants{

	const int nfreq = 1024;
	const int chunk_size = 1024 * 384 * 2; // 1024 intensiy samples per chunk, 384 baseband samples per intensity sample, 2 polarization
	const int max_chunks = 8;
	//const int max_processors = 10;
	const int frame_per_second = 390625;
	const unsigned int buffer_size = 819200; 
	const int file_packets = 131072; // vdif pakcets per vdif file
	const int udp_packets = 32; //change this to 8 on boris
	const int nfft = 16;
	const int nsample_integrate = 384; // 384 * 2.56us = 0.98 ms
	const int nt = chunk_size / nsample_integrate; // number of intensity samples per chunk
	int intensity_chunk_size = nt * nfreq;
	int intensity_buffer_size = intensity_chunk_size * max_chunks;
}


assembled_chunk::assembled_chunk() : t0(0) {

	data = new unsigned char[constants::chunk_size * (constants::nfreq+32)]; //additional 32 comes from header
}

void assembled_chunk::set_time(long unsigned int time) {
	t0 = time;
}

void assembled_chunk::set_data(int i, unsigned char x) {
	data[i] = x;
}


assembled_chunk::~assembled_chunk() {
	delete[] data;
}

//processor class
vdif_processor::vdif_processor(int id, bool flag){

	p_id = id;
	is_running = false;
	upch = flag;
	processor_chunk = new assembled_chunk();
	if (upch) {
		temp_buf = new int[(constants::nfreq+32)*constants::nfft];
		intensity_buf = new int[(constants::nfreq+32)*constants::nfft];
	} else {
		temp_buf = new int[constants::nfreq+32];
		intensity_buf = new int[constants::nfreq+32];
	}

	fft_buf = new float[constants::nsample_integrate * constants::nfreq];

	CPU_ZERO(&p_cpuset);
	CPU_SET((p_id % 4) + 4, &p_cpuset); //assign core 5,6,7,8 to do the square&sum
}

vdif_processor::~vdif_processor(){
	delete[] fft_buf;
	delete[] temp_buf;
	delete[] intensity_buf;
	delete processor_chunk;
}	

//processing method
void vdif_processor::process_chunk(int *intensity, int index, char &mask) {
	
	sched_setaffinity(0, sizeof(cpu_set_t), &p_cpuset);	

	int seconds,frames;
	long int t0;	
	
	memcpy(&seconds,processor_chunk->data,sizeof(int));
        memcpy(&frames,(processor_chunk->data) + 4,sizeof(int));
        t0 = (long int) (seconds & 0x3FFFFFFF) * (long int) constants::frame_per_second + (long int) (frames & 0xFFFFFF);

	cout << t0 << endl; //t0 defines a unique timestamp, can be used to check data stream continuity
	processor_chunk->t0 = t0;
	
	//mask-out invalid data
	for (int i = 0; i < constants::chunk_size; i++) {
		if (processor_chunk->data[i * 1056 + 3] >> 7) {
			for (int j = 0; j < 1056; j++) {
				processor_chunk->data[i*1056 + j] = 136; // 136 = 0+0j
			}
		}
	}

	for (int i = 0; i < constants::nt; i++) {
		//no upchannelization for now
		if (upch) {
			/*
			upchannelize_square_sum(constants::nsample_integrate, constants::nfreq + 32, constants::nfft, processor_chunk->data + i * constants::nsample_integrate * (constants::nfreq+32),fft_buf, p);
			for (int j = 0; j < constants::nfreq * constants::nfft; j++) {
				intensity_buf[j + 32 * constants::nfft] = 0;
				for (int k = 0; k < constants::nsample_integrate / constants::nfft; k++) {
					intensity_buf[j + 32 * constants::nfft] += (int) fft_buf[(j + 32 * constants::nfft) + (k * constants::nfreq + 32)  * constants::nfft];
				}
				intensity[j * constants::nt + i] = intensity_buf[j + 32 * constants::nfft];
			}
		*/	
		} else {
			//square & sum written by Andre
			//The '32' trick: treat header bytes as data, square & sum everything, and only took the last 1024 results
			u4_square_and_sum(constants::nsample_integrate, constants::nfreq + 32, processor_chunk->data + i * constants::nsample_integrate * (constants::nfreq+32), temp_buf, intensity_buf);
			//put results into intensity buffer
			for (int j = 0; j < constants::nfreq; j++) {
				intensity[j * constants::nt + i] = intensity_buf[j + 32];
			}
		}
	}


	cout << "Processing chunk done." << endl;
	
	//tell the assembler that I am done
 	unique_lock<mutex> lk3(mtx3);
	mask = mask | (1 << index);
	cv[p_id].notify_one();
	lk3.unlock();

	//ready for next job
	unique_lock<mutex> lk2(mtx2);
	is_running = false;
	lk2.unlock();

}


//assembler class
vdif_assembler::vdif_assembler(const char *arg1, const char *arg2, bool flag1, bool flag2, int n){
	
	number_of_processors = n;
	if (strcmp("network",arg1)==0) {	
		mode = 0;
		port = atoi(arg2);
	} else if (strcmp("disk",arg1)==0) {
		mode = 1;
		filelist_name = new char[strlen(arg2)];
		strcpy(filelist_name,arg2);
	} else if (strcmp("simulate",arg1)==0) {
		mode = 2;
	} else {
		cout << "Unsupported option." << endl;
		exit(1);
	}
	
	number_of_processors = n;
	processors = new vdif_processor *[number_of_processors];
	
	upch = flag1;
	for (int i = 0; i< number_of_processors; i++) {
		processors[i] = new vdif_processor(i, flag1);
	}

	if (upch) {
		fftw_complex *in, *out;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * constants::nfft);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * constants::nfft);
		p = fftw_plan_dft_1d(constants::nfft, in, out, FFTW_FORWARD, FFTW_MEASURE);
		fftw_free(in);fftw_free(out);
	}


	bufsize = 0; p_index = 0;

	chunk_count = 0;
	i_start_index = 0; i_end_index = 0;
	processor_threads = new thread[number_of_processors];
	cv = new condition_variable[number_of_processors];

	if (upch) {
		constants::intensity_chunk_size *= 16;
		constants::intensity_buffer_size *= 16; 
	}
	intensity_buffer = new int[constants::intensity_buffer_size];
	job_list = new int[number_of_processors];

    output = fopen("test.dat","w");
    if (!output) {
		cout << "Can't write to file." << endl;
		exit(10);
    }
	write_to_disk = flag2;


}

vdif_assembler::~vdif_assembler() {
	
	fclose(output);

	delete[] processors;
}

void vdif_assembler::run() {

	thread stream_t;

	if (mode==0) {
		stream_t = thread(&vdif_assembler::network_capture,this);
	} else if (mode==1) {
		stream_t = thread(&vdif_assembler::read_from_disk,this);
	} else if (mode==2) {
		stream_t = thread(&vdif_assembler::simulate,this);
	}
	
	stream_t.join();

}

//this method feeds bonsai an intensity chunk
void vdif_assembler::get_intensity_chunk(float *intensity, ssize_t stride) {

	//check if the chunk is ready
	unique_lock<mutex> lk3(mtx3);		
	if (!(intensity_buffer_mask & (1 << i_start_index))) {
		cv[job_list[i_start_index]].wait(lk3);
	}
		
	//copy data over
	for (int i = 0; i < constants::nt; i++ ) {
		for (int j = 0; j < constants::nfreq; j++) {
		
		intensity[i*stride + j] = intensity_buffer[i_start_index * constants::intensity_chunk_size + i * constants::nfreq + j];
				
		}
	}
	chunk_count++;

	//write 32-bits data to disk,can be changed to 8-bit
	if (write_to_disk) {	
		fwrite(intensity_buffer + i_start_index * constants::intensity_chunk_size, sizeof(int), constants::intensity_chunk_size, output);
		cout << "Writing chunk " << chunk_count << " to test.dat..." << endl;
	}

	intensity_buffer_mask = intensity_buffer_mask ^ (1 << i_start_index);	
	lk3.unlock();
		
	//measure performance
	perf = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count();
	cout << perf << "ms" <<endl;
	start = chrono::steady_clock::now();
	tot_perf += perf;
	avg_perf = tot_perf / chunk_count;
	cout << "Average performance: " << avg_perf << "ms" << endl;

	cout << "Finish assembling chunk " << chunk_count << ", sending to bonsai..." << endl;
	i_start_index = (i_start_index + 1) % constants::max_chunks;

}

//job distributor
void vdif_assembler::assign_chunk() {

	cout << "Assigned to processor: " << p_index << endl;
    processor_threads[p_index] = thread(&vdif_processor::process_chunk,processors[p_index],intensity_buffer+i_end_index*constants::intensity_chunk_size,i_end_index,ref(intensity_buffer_mask));
    processors[p_index]->is_running = true;
    processor_threads[p_index].detach();
	job_list[i_end_index] = p_index;
    i_end_index = (i_end_index + 1) % constants::max_chunks;
    bufsize = 0;
    p_index = get_free_processor();


}

//free-worker searcher
int vdif_assembler::get_free_processor() {
	int i = 0;

	//this is a bad busy loop, but not too bad since there will always be a free worker
	for (;;) {
		unique_lock<mutex> lk2(mtx2);
		if (!processors[i]->is_running) {
			return i;
		}
		i++;
		if (i == number_of_processors) {
			//cout << "All processors are busy. Waiting for them to finish..." << endl;
			i = 0;
		}
		lk2.unlock();
	}

}

void vdif_assembler::network_capture() {
	
	cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(5,&cpuset); //assign core 6 to this thread since 6 is a perfect number
    sched_setaffinity(0,sizeof(cpu_set_t),&cpuset);

	
	int size = constants::udp_packets * 1056; //every packet has 1056 bytes
	struct sockaddr_in server_address;
	memset(&server_address, 0, sizeof(server_address));
	
	int sock_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (sock_fd < 0) {
		cout << "socket failed." << std::endl;
	}
	server_address.sin_family = AF_INET;
	server_address.sin_port = htons(port);
	server_address.sin_addr.s_addr = inet_addr("127.0.0.1"); //change this to 10.1.1.2 on boris
	if (bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address)) < 0) {
		cout << "bind failed." << std::endl;
	}
	
	for (;;) {
		if (read(sock_fd, (processors[p_index]->processor_chunk->data) + bufsize * 1056, size) == size) {
			
			bufsize += constants::udp_packets;
			if (bufsize == constants::chunk_size) {
				//cout << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count() << "ms" <<endl;
		                //start = chrono::steady_clock::now();
				this->assign_chunk();

			}
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
			fread(processors[p_index]->processor_chunk->data + bufsize * 1056,sizeof(char),1056,fp);
			bytes_read += 1056;
			bufsize++;
			if (bufsize == constants::chunk_size) {
				this->assign_chunk();
			}

		}
		
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
						processors[p_index]->processor_chunk->data[bufsize*1056 + i*4+j] =(int) ((word[i] >> (j*8)) & 0xFF);
					}
				}

					
				for (int i = 32; i < 1056; i++) {
					processors[p_index]->processor_chunk->data[bufsize*1056+i] = voltage;
				}
	
				bufsize ++;
				if (bufsize == constants::chunk_size) {
					this->assign_chunk();
				}
                        	
				
			}	
		}
		//this_thread::sleep_for(chrono::seconds(1));
	}

}

	
