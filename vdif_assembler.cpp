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
#include <ipp.h>
#include <rf_pipelines.hpp>
#include <sys/time.h>
#include <chrono>
//#include "gaussian.cpp"

#include <sched.h>
#include <pthread.h>

#include "vdif_assembler.hpp"
#include "square_sum_andre.cpp"
//#include "upchannelize.cpp"
#include <assert.h>


using namespace std;
using namespace rf_pipelines;

auto start = chrono::steady_clock::now();
mutex mtx,mtx2,mtx3;
condition_variable cv,cv_chunk;

float perf, tot_perf, avg_perf;
//fftwf_plan p;


char *strptime(const char * __restrict, const char * __restrict, struct tm * __restrict);

namespace constants{

	const int nfreq = 1024;
	const int header_size = 12; //int64 time + int32 thread ID
	const int chunk_size = 1024 * 384 * 2;
	const int max_chunks = 8;
	//const int max_processors = 10;
	const int frame_per_second = 390625;
	const unsigned int buffer_size = 819200; //
	const int file_packets = 131072;
	const int udp_packets = 32;
	const int nfft = 16;
	const int nsample_integrate = 16*24*2;
	const int nt = chunk_size / nsample_integrate;
	const int intensity_chunk_size = nt * nfreq;
	const int intensity_buffer_size = intensity_chunk_size * max_chunks;
}


assembled_chunk::assembled_chunk() : t0(0) {

	data = new unsigned char[constants::chunk_size * (constants::nfreq+32)];

}

void assembled_chunk::set_time(long unsigned int time) {
	t0 = time;
}

void assembled_chunk::set_data(int i, unsigned char x) {
	//assert(i < constants::chunk_size * constants::nfreq);
	//assert(i >= 0);
	data[i] = x;
}


assembled_chunk::~assembled_chunk() {
	delete[] data;
}

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

	CPU_ZERO(&p_cpuset);
	CPU_SET(p_id * 2 + 3, &p_cpuset);
}

vdif_processor::~vdif_processor(){
	delete[] temp_buf;
	delete[] intensity_buf;
	delete processor_chunk;
}

void vdif_processor::process_chunk(int *intensity, int index, char &mask) {
	
	sched_setaffinity(0, sizeof(cpu_set_t), &p_cpuset);	

	int seconds,frames;
	long int t0;	
	
	memcpy(&seconds,processor_chunk->data,sizeof(int));
        memcpy(&frames,(processor_chunk->data) + 4,sizeof(int));
        t0 = (long int) (seconds & 0x3FFFFFFF) * (long int) constants::frame_per_second + (long int) (frames & 0xFFFFFF);

	cout << t0 << endl;
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
		if (upch) {
			
			//upchannelize(constants::nsample_integrate, constants::nfreq + 32, constants::nfft, processor_chunk->data + i * constants::nsample_integrate * (constants::nfreq+32));
			
		} else {
			u4_square_and_sum(constants::nsample_integrate, constants::nfreq + 32, processor_chunk->data + i * constants::nsample_integrate * (constants::nfreq+32), temp_buf, intensity_buf);
			for (int j = 0; j < constants::nfreq; j++) {
				intensity[j * constants::nt + i] = intensity_buf[j + 32];
			}
		}
	}

	//gaussian_number(intensity, 1024*1024);

	cout << "Processing chunk done." << endl;

	unique_lock<mutex> lk3(mtx3);
 
	mask = mask | (1 << index);
	
	lk3.unlock();

	unique_lock<mutex> lk2(mtx2);
	is_running = false;
	lk2.unlock();

}


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

	bufsize = 0; p_index = 0;

	chunk_count = 0;
	i_start_index = 0; i_end_index = 0;
	processor_threads = new thread[number_of_processors];
	intensity_buffer = new int[constants::intensity_buffer_size];
	intensity_buffer_mask = 0;

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

/*
int vdif_assembler::is_full() {

	return bufsize == constants::buffer_size;
}
*/

void vdif_assembler::run() {

	//thread assemble_t(&vdif_assembler::assemble_chunk,this);
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
	//assemble_t.join();

}
/*
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
*/
void vdif_assembler::get_intensity_chunk(float *intensity, ssize_t stride) {

	for (;;) {
		unique_lock<mutex> lk3(mtx3);		
		if (intensity_buffer_mask & (1 << i_start_index)) {

			for (int i = 0; i < constants::nt; i++ ) {
				for (int j = 0; j < constants::nfreq; j++) {
		
				intensity[i*stride + j] = intensity_buffer[i_start_index * constants::intensity_chunk_size + i * constants::nfreq + j];
				
				}
			}
			chunk_count++;
			if (write_to_disk) {	
				fwrite(intensity_buffer + i_start_index * constants::intensity_chunk_size, sizeof(int), constants::intensity_chunk_size, output);
                        	cout << "Writing chunk " << chunk_count << " to test.dat..." << endl;
			}
			
			intensity_buffer_mask = intensity_buffer_mask ^ (1 << i_start_index);
			perf = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count();
			cout << perf << "ms" <<endl;
			start = chrono::steady_clock::now();
			tot_perf += perf;
			avg_perf = tot_perf / chunk_count;
			cout << "Average performance: " << avg_perf << "ms" << endl;
			lk3.unlock();
			cout << "Finish assembling chunk " << chunk_count << ", sending to bonsai..." << endl;
			i_start_index = (i_start_index + 1) % constants::max_chunks;
			break;
		}
		lk3.unlock();
	}



}
/*
void vdif_assembler::assemble_chunk() {

	int p_index;
	bool assigned;
	long int t0;
	int seconds;
	int frames;
	for (;;) {

		//cout << " start: " << start_index << " end: " << end_index << endl;
		unique_lock<mutex> lk(mtx);
		
		if (bufsize < constants::chunk_size) {
			cv.wait(lk);
		}
		bufsize -= constants::chunk_size;
		cout << "Chunk found. excess: " << bufsize << endl;
		lk.unlock();
		cout << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count() << "ms" <<endl;
                start = chrono::steady_clock::now();
		assigned = false;
		p_index = 0;
		memcpy(&seconds,data_buf,sizeof(int));
		memcpy(&frames,data_buf+4,sizeof(int));
		t0 = (long int) (seconds & 0x3FFFFFFF) * (long int) constants::frame_per_second + (long int) (frames & 0xFFFFFF);
		
		while (!assigned) { 
			unique_lock<mutex> lk2(mtx2);
			if (!processors[p_index]->is_running) {
				
				processors[p_index]->processor_chunk->set_time(t0);
				for (int i = 0; i < constants::chunk_size; i++) {
					for (int j = 0; j < constants::nfreq; j++) {
						processors[p_index]->processor_chunk->set_data(i*constants::nfreq + j,data_buf[(start_index+i)*constants::nfreq + i * 4 + j]);
					}
				}
				start_index = (start_index + constants::chunk_size) % constants::buffer_size;
				
				cout << "Assigned to processor: " << p_index << endl;
 	
				processor_threads[p_index] = thread(&vdif_processor::process_chunk,processors[p_index],intensity_buffer+i_end_index*constants::intensity_chunk_size,i_end_index,ref(intensity_buffer_mask));
				processor_threads[p_index].detach();
				assigned = true;
				i_end_index = (i_end_index + 1) % constants::max_chunks;
				
			} else {
				p_index++;
			}
			lk2.unlock();
			if (p_index == number_of_processors) {
				cout << "All processors are busy. Waiting for them to finish..." << endl;
				p_index = 0;
			}

			
		}
				
	}
}
*/

void vdif_assembler::assign_chunk() {

	cout << "Assigned to processor: " << p_index << endl;

        processor_threads[p_index] = thread(&vdif_processor::process_chunk,processors[p_index],intensity_buffer+i_end_index*constants::intensity_chunk_size,i_end_index,ref(intensity_buffer_mask));
        processors[p_index]->is_running = true;
        processor_threads[p_index].detach();

        i_end_index = (i_end_index + 1) % constants::max_chunks;
        bufsize = 0;
        p_index = get_free_processor();


}
int vdif_assembler::get_free_processor() {
	int i = 0;
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
        CPU_SET(6,&cpuset);
        sched_setaffinity(0,sizeof(cpu_set_t),&cpuset);

	
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


/*
void vdif_assembler::vdif_read(unsigned char *data, int size) {


	int word[8];
	int count = 0;
	long int t0 = 0;
	int pol = 0;
	int nmissing = 0;
	long int current,expect;
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
		int prev = (end_index - 1) % constants::buffer_size;
		expect = header_buf[prev].t0 * 2 + header_buf[prev].polarization + 1;
		//cout << current << " " << expect << endl;
		if (expect!=1) {
			//nmissing = (int) (current - expect);
		}

		if (nmissing) {
			
			cout << "Missing "<< nmissing << " packets." << endl;
			cout << "current: " << current << " expect: " << expect << endl; 
			cout << "start: " << start_index << " end: " << end_index << endl;
			fill_missing(nmissing);
		}

		//expect = current + 1;
		
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
*/

	
