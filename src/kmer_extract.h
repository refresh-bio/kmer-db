#pragma once
#include "alphabet.h"

#include <vector>
#include <memory>
#include <cinttypes>

class KmerHelper {

public:

	template <class Filter>
	static size_t extract(
		char* sequence,
		size_t sequenceLength,
		uint32_t kmerLength,
		const Alphabet& alphabet,
		const Filter& filter,
		kmer_t* kmers) {

	
		size_t counter = 0;

		kmer_t kmer_str, kmer_rev, kmer_can;
		uint32_t kmer_len_shift = (kmerLength - 1) * alphabet.bitsPerSymbol;
		kmer_t kmer_mask = (1ull << (alphabet.bitsPerSymbol * kmerLength)) - 1;
		int omit_next_n_kmers;
		uint32_t i;

		kmer_str = kmer_rev = 0;

		uint32_t str_pos = kmer_len_shift - alphabet.bitsPerSymbol;
		uint32_t rev_pos = alphabet.bitsPerSymbol;

		omit_next_n_kmers = 0;

		// calculate k-mers shifting to get prefix of at least 8 bits
		size_t kmer_prefix_shift = 0;
		kmer_t tail_mask = 0;
		int prefix_bits = (int)kmerLength * alphabet.bitsPerSymbol - SUFFIX_BITS;

		if (prefix_bits < 8) {
			kmer_prefix_shift = (size_t)(8 - prefix_bits);
			tail_mask = (1ULL << kmer_prefix_shift) - 1;
		}


		for (i = 0; i < kmerLength - 1; ++i, str_pos -= alphabet.bitsPerSymbol, rev_pos += alphabet.bitsPerSymbol)
		{
			int8_t symb = alphabet.map(sequence[i]);
			if (symb < 0)
			{
				symb = 0;
				omit_next_n_kmers = i + 1;
			}
			kmer_str += (kmer_t)symb << str_pos;
			kmer_rev += (kmer_t)(3 - symb) << rev_pos; // this makes sense only for DNA alphabet
		}

		for (; i < sequenceLength; ++i)
		{
			int8_t symb = alphabet.map(sequence[i]);
			if (symb < 0)
			{
				symb = 0;
				omit_next_n_kmers = kmerLength;
			}
			kmer_str = (kmer_str << alphabet.bitsPerSymbol) + (kmer_t)symb;
			kmer_str &= kmer_mask;

			kmer_rev >>= alphabet.bitsPerSymbol;
			kmer_rev += (kmer_t)(alphabet.size - 1 - symb) << kmer_len_shift;

			if (omit_next_n_kmers > 0)
			{
				--omit_next_n_kmers;
				continue;
			}

			if (alphabet.preserveStrand) {
				kmer_can = kmer_str;
			}
			else {
				kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;
			}

			// ensure at least 8-bit prefix
			kmer_can = (kmer_can << kmer_prefix_shift) | (kmer_can & tail_mask);

			if (filter(kmer_can)) {
				kmers[counter++] = kmer_can;
			}

		}

		return counter;
	}

	static void sort(kmer_t* kmers, size_t count, uint32_t n_threads = 1) {
//		ParallelSort(kmers, count, n_threads);
		refresh::sort::pdqsort_branchless(kmers, kmers + count);
	}

	static void sortAndUnique(kmer_t* kmers, size_t& count, uint32_t n_threads = 1) {
//		ParallelSort(kmers, count, n_threads);
		refresh::sort::pdqsort_branchless(kmers, kmers + count);
		auto it = std::unique(kmers, kmers + count);

		count = it - kmers;
	}

	static void unique(kmer_t* kmers, size_t& count) {
//		std::sort(kmers, kmers + count);
		refresh::sort::pdqsort_branchless(kmers, kmers + count);
		auto it = std::unique(kmers, kmers + count);

		count = it - kmers;
	}

};

