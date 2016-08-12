/// @brief Squares and sums 4 bit data with the form
/// [f0t0, f1t0, ... fNt0, f0t1, f1t1, ... fNt1, ... f0tM, f1TM, ... fNtM]
/// f is frequency, t is time, N is number of frequencies, M is number of samples to integration.
/// Each value is a (4+4)-bit offset encoded complex number.
/// @param intergration_time The number of time samples to sum (M)
/// @param num_freq The number of frequencies in each time sample (N)
/// @parma data The blob of data with the format above.
/// @parap temp_buf A working array of size = num_freq * sizeof(int)
/// @param sum The results of the square and sum; size = num_freq * sizeof(int)
/// @returns The sum array contains the square and sum of the values in the order [f0, f1, f2, бн, fN]
#include <immintrin.h>

inline void u4_square_and_sum(const int integration_time, const int num_freq,
	const unsigned char * data, int * temp_buf, int * sum) {

	for (int frame = 0; frame < integration_time; frame++) {
		for (int freq = 0; freq < num_freq / 32; ++freq) {

			const int index = frame * num_freq + freq * 32;

			__m256i ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;

			// Load 64 4 bit numbers
			ymm0 = _mm256_loadu_si256((__m256i const *)&data[index]);

			// Shift the high 4-bits to the low 4-bits in each 8 bit block
			ymm1 = _mm256_srli_epi64(ymm0, 4); // real

											   // Mask out the lower 4 bits
			ymm2 = _mm256_set1_epi32(0x0f0f0f0f);
			ymm0 = _mm256_and_si256(ymm0, ymm2); // imag
			ymm1 = _mm256_and_si256(ymm1, ymm2); // real

												 // This packs (real, imag) (8+8) pairs together
			ymm3 = _mm256_unpacklo_epi8(ymm0, ymm1);
			ymm4 = _mm256_unpackhi_epi8(ymm0, ymm1);

			// subtract 8 to make the 8-bit numbers twos complement
			ymm2 = _mm256_set1_epi8(8);
			ymm3 = _mm256_sub_epi8(ymm3, ymm2);
			ymm4 = _mm256_sub_epi8(ymm4, ymm2);

			// Take the abs value since the multi is unsigned
			ymm3 = _mm256_abs_epi8(ymm3);
			ymm4 = _mm256_abs_epi8(ymm4);

			// Multiply and add real and imaginary 8+8-bit pairs into 16-bit ints
			ymm5 = _mm256_maddubs_epi16(ymm3, ymm3); // hi
			ymm6 = _mm256_maddubs_epi16(ymm4, ymm4); // lo

													 // Extend to 32-bit
			ymm7 = _mm256_set1_epi32(0);
			ymm0 = _mm256_unpacklo_epi16(ymm5, ymm7);
			ymm1 = _mm256_unpackhi_epi16(ymm5, ymm7);
			ymm2 = _mm256_unpacklo_epi16(ymm6, ymm7);
			ymm3 = _mm256_unpackhi_epi16(ymm6, ymm7);

			int out_index = freq * 32;

			if (frame != 0) {
				ymm4 = _mm256_loadu_si256((__m256i const *)&temp_buf[out_index + 0 * 8]);
				ymm5 = _mm256_loadu_si256((__m256i const *)&temp_buf[out_index + 1 * 8]);
				ymm6 = _mm256_loadu_si256((__m256i const *)&temp_buf[out_index + 2 * 8]);
				ymm7 = _mm256_loadu_si256((__m256i const *)&temp_buf[out_index + 3 * 8]);

				ymm0 = _mm256_add_epi32(ymm0, ymm4);
				ymm1 = _mm256_add_epi32(ymm1, ymm5);
				ymm2 = _mm256_add_epi32(ymm2, ymm6);
				ymm3 = _mm256_add_epi32(ymm3, ymm7);
			}

			_mm256_storeu_si256((__m256i *)&temp_buf[out_index + 0 * 8], ymm0);
			_mm256_storeu_si256((__m256i *)&temp_buf[out_index + 1 * 8], ymm1);
			_mm256_storeu_si256((__m256i *)&temp_buf[out_index + 2 * 8], ymm2);
			_mm256_storeu_si256((__m256i *)&temp_buf[out_index + 3 * 8], ymm3);
		}
	}

	// Reorder the numbers.
	for (int i = 0; i < num_freq; ++i) {
		// Fix stupid index problem
		int m32 = i % 32;
		if (m32 < 16) {
			m32 = (m32 / 4) * 4;
		}
		else {
			m32 = -12 + ((m32 - 16) / 4) * 4;
		}
		// TODO this should be a non-temporal store.
		sum[i] = temp_buf[i + m32];
	}
}