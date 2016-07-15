#include "vdif_assembler.cpp"
using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 4) {
		cout << "run with:test <network/disk/simulate> <port/filelist> <-u/-nu>" << endl;
		exit(10);
	}
	int n = 12;
	bool up = false;
	if (strcmp("-u",argv[3])==0) {
		cout << "Upchannelization enabled." << endl;
		up = true;
	}
	vdif_assembler a(argv[1],argv[2],up,n);
	//vdif_processor p[n];
	
	a.run();
}
