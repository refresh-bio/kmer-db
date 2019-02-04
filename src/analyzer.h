#pragma once
#include "kmer_db.h"
#include "filter.h"


class Analyzer {
public:

	std::unique_ptr<SetFilter> selectSeedKmers(const FastKmerDb& db, size_t minSamplesCount);

	void printStats(const FastKmerDb& db);

	void operator()(const FastKmerDb& db);
};