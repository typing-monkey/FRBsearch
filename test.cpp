#include "vdif_assembler.cpp"
using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cout << "run with:test <network/disk/simulate> <port/filelist> " << endl;
		exit(10);
	}

	vdif_assembler a(argv[2],argv[3]);
	vdif_processor p(true, "test.dat");
	a.register_processor(&p);
	a.run();
}

