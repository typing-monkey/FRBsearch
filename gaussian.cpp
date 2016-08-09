#include <random>

using namespace std;

void gaussian_number(int *data, int n) {

        double sample_rms = 1.0;
        random_device rd;
        mt19937 rng(rd());
        normal_distribution<float> dist(0, sample_rms);

        for (int i = 0; i < n; i++) {
                data[i] = (int) dist(rng);
        }


}

