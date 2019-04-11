#pragma once
#include "prefix_kmer_db.h"
#include "filter.h"


class Analyzer {
public:

	std::unique_ptr<SetFilter> selectSeedKmers(const PrefixKmerDb& db, size_t minSamplesCount);

	void printStats(const PrefixKmerDb& db);

	void operator()(const PrefixKmerDb& db);
};