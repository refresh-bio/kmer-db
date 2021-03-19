#pragma once
#include <vector>
#include <memory>

template <class Filter>
size_t extractKmers(
	char* sequence,
	size_t sequenceLength,
	uint32_t kmerLength,
	std::shared_ptr<Filter> filter,
	kmer_t* kmers,
	uint32_t* positions) {

	
	static char* map = []() {
		static char _map[256];
		std::fill_n(_map, 256, -1);
		_map['a'] = _map['A'] = 0;
		_map['c'] = _map['C'] = 1;
		_map['g'] = _map['G'] = 2;
		_map['t'] = _map['T'] = 3;
		return _map;
	}();
	
	size_t counter = 0;

	kmer_t kmer_str, kmer_rev, kmer_can;
	uint32_t kmer_len_shift = (kmerLength - 1) * 2;
	kmer_t kmer_mask = (1ull << (2 * kmerLength)) - 1;
	int omit_next_n_kmers;
	uint32_t i;
	
	kmer_str = kmer_rev = 0;

	uint32_t str_pos = kmer_len_shift - 2;
	uint32_t rev_pos = 2;

	omit_next_n_kmers = 0;

	// calculate k-mers shifting to get prefix of at least 8 bits
	size_t kmer_prefix_shift = 0;
	kmer_t tail_mask = 0;
	int prefix_bits = ((int)kmerLength - SUFFIX_LEN) * 2;
	
	if (prefix_bits < 8) {
		kmer_prefix_shift = (size_t)(8 - prefix_bits);
		tail_mask = (1ULL << kmer_prefix_shift) - 1;
	}


	for (i = 0; i < kmerLength - 1; ++i, str_pos -= 2, rev_pos += 2)
	{
		char symb = map[static_cast<unsigned char>(sequence[i])];
		if (symb < 0)
		{
			symb = 0;
			omit_next_n_kmers = i + 1;
		}
		kmer_str += (kmer_t)symb << str_pos;
		kmer_rev += (kmer_t)(3 - symb) << rev_pos;
	}

	for (; i < sequenceLength; ++i)
	{
		char symb = map[static_cast<unsigned char>(sequence[i])];
		if (symb < 0)
		{
			symb = 0;
			omit_next_n_kmers = kmerLength;
		}
		kmer_str = (kmer_str << 2) + (kmer_t)symb;
		kmer_str &= kmer_mask;

		kmer_rev >>= 2;
		kmer_rev += (kmer_t)(3 - symb) << kmer_len_shift;

		if (omit_next_n_kmers > 0)
		{
			--omit_next_n_kmers;
			continue;
		}

		kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;

		// ensure at least 8-bit prefix
		kmer_can = (kmer_can << kmer_prefix_shift) | (kmer_can & tail_mask);

	//	filter->add(kmer_can);	
		if ((*filter)(kmer_can)) {
			
			if (positions != nullptr) {
				positions[counter] = i;
			}
			
			kmers[counter++] = kmer_can;
		}
		
	}

	return counter;
}