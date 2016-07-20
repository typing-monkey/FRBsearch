#include "aro_stream.cpp"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 4) {
		cout << "run with:test <network/disk/simulate> <port/filelist> <-u/-nu>" << endl;
		exit(10);
	}

	int n = 12;
	bool up = false;
	
	int nfreq = 1024;
	double freq_lo_MHz = 400.0;
	double freq_hi_MHz = 800.0;
	double dt_sample = 6.5536e-4;

	if (strcmp("-u",argv[3])==0) {
		cout << "Upchannelization enabled." << endl;
		up = true;
	}
	
	shared_ptr<wi_stream> stream = make_shared<aro_stream>(nfreq,freq_lo_MHz,freq_hi_MHz,dt_sample, argv[1], argv[2], up, n);

	vector<shared_ptr<wi_transform>> transform_list;
	
	string bonsai_config_filename = "bonsai_config.hdf5";
	string bonsai_output_filename = "bonsai_outputs.hdf5";	

	transform_list.push_back(make_bonsai_dedisperser(bonsai_config_filename,bonsai_output_filename));
	stream->run(transform_list);

	//vdif_assembler a(argv[1],argv[2],up,n);
	//vdif_processor p[n];
	
	//a.run();
}

