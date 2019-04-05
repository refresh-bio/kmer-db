#pragma once
#include "prefix_kmer_db.h"

class SimilarityCalculator {
public:
	SimilarityCalculator(int _num_threads, size_t cacheBufferMb);

	virtual void operator()(PrefixKmerDb& db, LowerTriangularMatrix<uint32_t>& matrix) const;

	virtual void operator()(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<uint32_t>& vector) const;


protected:

	static const int PREFETCH_DIST = 48;

	int num_threads;

	size_t cacheBufferMb;

	bool avx2_present;
};