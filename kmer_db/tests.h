#pragma once
#include "kmer_db.h"


class Tests {
public:
	static void compareWithNaive(const string& fname);


	static void comparePatterns(const AbstractKmerDb& db1, const AbstractKmerDb& db2, const std::vector<kmer_t>& kmers);
	static void testSerialization(const FastKmerDb& db);
	static void testDistanceMatrix(const AbstractKmerDb& db, const std::string& fname);
};