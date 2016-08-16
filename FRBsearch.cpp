#include "aro_stream.cpp"
#include "preprocessing.cpp"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 5) {
		cout << "run with:FRBsearch <network/disk/simulate> <port/filelist> <-u/-nu> <-w/-nw>" << endl;
		exit(10);
	}
	
	
	int n = 3;
	bool up = false;
	bool write = false;
	int nfreq = 1024;
	double freq_lo_MHz = 400.0;
	double freq_hi_MHz = 800.0;
	double dt_sample = 9.8304e-4;

	if (strcmp("-u",argv[3])==0) {
		cout << "Upchannelization is not implemented yet." << endl;
		exit(10);
		cout << "Upchannelization enabled." << endl;
		up = true;
	}
	
	if (strcmp("-w",argv[4])==0) {
		cout << "Writing to intensity data to test.dat." << endl;
		write = true;
	}

	shared_ptr<wi_stream> stream = make_shared<aro_stream>(nfreq,freq_lo_MHz,freq_hi_MHz,dt_sample, argv[1], argv[2], up, write, n);

	vector<shared_ptr<wi_transform>> transform_list;
	
	string bonsai_config_filename = "bonsai_config.hdf5";
	string bonsai_output_filename = "bonsai_outputs.hdf5";	
	int nt_per_file = 16384;

	//rfi filters

	shared_ptr<wi_transform> stb = shared_ptr<wi_transform>(new sys_temperature_bandpass());
	shared_ptr<wi_transform> r_o1 = shared_ptr<wi_transform>(new remove_outliers(5.0));
	//insert hp_filter here

	double hpf_width = 0.2;
	shared_ptr<wi_transform> hpf = shared_ptr<wi_transform>(new highpass_filter((int) hpf_width/dt_sample));
	shared_ptr<wi_transform> r_o2 = shared_ptr<wi_transform>(new remove_outliers(5.0));
	shared_ptr<wi_transform> r_nf1 = shared_ptr<wi_transform>(new remove_noisy_freq(3.0));
	shared_ptr<wi_transform> r_bt = shared_ptr<wi_transform>(new remove_bad_times(2.0));
	shared_ptr<wi_transform> r_c = shared_ptr<wi_transform>(new remove_continuum());
	shared_ptr<wi_transform> r_nf2 = shared_ptr<wi_transform>(new remove_noisy_freq(3.0));

	/*
	transform_list.push_back(stb);
	transform_list.push_back(r_o1);
	transform_list.push_back(hpf);
	transform_list.push_back(r_o2);
	transform_list.push_back(r_nf1);
	transform_list.push_back(r_bt);
	transform_list.push_back(r_c);
	transform_list.push_back(r_nf2);
	*/
	transform_list.push_back(make_bonsai_dedisperser(bonsai_config_filename,bonsai_output_filename,nt_per_file));

	stream->run(transform_list);

	//vdif_assembler a(argv[1],argv[2],up,n);
	//vdif_processor p[n];
	
	//a.run();
}

