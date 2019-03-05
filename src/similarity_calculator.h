#pragma once
#include "kmer_db.h"

class SimilarityCalculator {
public:
	SimilarityCalculator(int _num_threads, size_t cacheBufferMb);

	virtual void operator()(FastKmerDb& db, LowerTriangularMatrix<uint32_t>& matrix) const;

	virtual void operator()(const FastKmerDb& db, const std::vector<kmer_t>& kmers, std::vector<uint32_t>& vector) const;


protected:
	int num_threads;

	size_t cacheBufferMb;

	bool avx2_present;
};