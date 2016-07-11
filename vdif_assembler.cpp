#include <time.h>
#include <sys/socket.h>
#include <fstream>
#include <stdio.h>
#include <queue>
#include <stdlib.h>
#include <errno.h>
#include <arpa/inet.h>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <fftw3.h>

#include "vdif_assembler.hpp"
#include "square_sum.cpp"
#include "upchannelize.cpp"

using namespace std;

mutex mtx;
condition_variable cv;

char *strptime(const char * __restrict, const char * __restrict, struct tm * __restrict);

namespace constants{

	const int nfreq = 1024;
	const int header_size = 12; //int64 time + int32 thread ID
	const int chunk_size = 98304;
	const int max_processors = 10;
	const int frame_per_second = 390625;
	const int buffer_size = 524288; //at most 8 chunk in buffer
	const int file_packets = 131072;
	const int udp_packets = 32;
	const int max_chunks = 10;
	const int nfft = 16;
}


assembled_chunk::assembled_chunk(long int start_time) {
	
	data = new unsigned char[constants::chunk_size * constants::nfreq];
	t0 = start_time;

}

assembled_chunk::~assembled_chunk() {
	delete[] data;
}

void assembled_chunk::set_data(int i, unsigned char x) {
	data[i] = x;
}

vdif_processor::vdif_processor(bool flag1, const char *filename){
	upch = flag1;
	output = fopen(filename, "w");
	if (!output) {
		cout << "Can't write to the file." << endl;
		exit(10);
	}	
	
}

vdif_processor::~vdif_processor(){
	fclose(output);
}

void vdif_processor::process_chunk(assembled_chunk *c) {
	cout << c->t0 << endl;	
	unsigned char temp[constants::chunk_size];

	for (int i = 0; i < constants::nfreq; i++) {
			for (int j = 0; j < (constants::chunk_size) >> 1; j++) {
				temp[j] = c->data[i + ((constants::nfreq*j) << 1)];
			}
			for (int j = 0; j < (constants::chunk_size) >> 1; j++) {
				temp[j+(constants::chunk_size>>1)] = c->data[i + ((constants::nfreq*j) << 1) + constants::nfreq];
			}
			int intensity[constants::chunk_size >> 8] = { 0 };
			int nintegrate = 8;
			if (upch) {
				upchannelize_square_sum(constants::chunk_size, constants::nfft, temp);
				nintegrate = 4;
			}
			square_sum(temp, constants::chunk_size, nintegrate, intensity);
			
			//for (int j = 0; j < int(constants::chunk_size / 256); j++) {
			//	if (intensity[j]) {
			//		cout << intensity[j] << endl;
			//	}
			//}
			fwrite(intensity, sizeof(intensity[0]), constants::chunk_size >> 8, output);
			//cout << "Frequency " << i << " : " << intensity[0] << endl;

		}
	
	cout << "Processing chunk done." << endl;

}


vdif_assembler::vdif_assembler(const char *arg1, const char *arg2){

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

	processors = new vdif_processor *[constants::max_processors];
	number_of_processors = 0;
	data_buf = new unsigned char[constants::buffer_size * constants::nfreq];
	header_buf = new struct header[constants::buffer_size];
	bufsize = 0;
	start_index = 0;
	end_index = 0;
	processor_threads = new thread[constants::max_processors];
	

}

vdif_assembler::~vdif_assembler() {
	delete[] processors;
	delete[] data_buf;
	delete[] header_buf;
}

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
}

int vdif_assembler::is_full() {

	return bufsize == constants::buffer_size;
}


void vdif_assembler::run() {

	thread assemble_t(&vdif_assembler::assemble_chunk,this);
	thread stream_t;

	if (mode==0) {
		stream_t = thread(&vdif_assembler::network_capture,this);
	} else if (mode==1) {
		stream_t = thread(&vdif_assembler::read_from_disk,this);
	} else if (mode==2) {
		stream_t = thread(&vdif_assembler::simulate,this);
	}
	
	stream_t.join();
		
	assemble_t.join();

}


void vdif_assembler::assemble_chunk() {
	
	for (;;) {
		//cout << " start: " << start_index << " end: " << end_index << endl;
		unique_lock<mutex> lk(mtx);
		
		if (bufsize < constants::chunk_size) {
			cv.wait(lk);
		}

		cout << "Chunk found" << endl;
		assembled_chunk c(header_buf[start_index].t0);	
		
		for (int i = 0; i < constants::chunk_size * constants::nfreq; i++) {
			c.set_data(i, data_buf[start_index*constants::nfreq+i]);
		}
		if (chunks.size() < constants::max_chunks) {
			chunks.push(&c);
		}
		start_index = (start_index + constants::chunk_size) % constants::buffer_size;

		bufsize -= constants::chunk_size;
	
		cout << "excess: " << bufsize << endl;
		
		for (int i = 0; i < number_of_processors; i++) {
			processors[i]->process_chunk(get_chunk());
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

assembled_chunk* vdif_assembler::get_chunk() {
	
	assembled_chunk* temp;
	if (chunks.size()) {
		temp = chunks.front();
		chunks.pop();
		return temp;
	} else {
		cout << "no chunk available at this moment." << endl;
		return 0;
	}

}

void vdif_assembler::simulate() {
	
	int word[8];
	struct tm epoch;
	strptime("2000-01-01 00:00:00","%Y-%m-%dT %H:%M:%S", &epoch);
	
	
	for(int t0 = difftime(time(0),mktime(&epoch));;t0++) {
		int duration = (int)((rand() % 3 + 3)*1000);
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
			nmissing = (int) (expect - current);
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


