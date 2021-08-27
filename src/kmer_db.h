#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "types.h"
#include "pattern.h"
#include "hashmap_lp.h"
#include "array.h"
#include "queue.h"
#include "aligned_vector.h"
#include "row_add.h"
#include "parallel_sorter.h"

#include <map>
#include <fstream>
#include <chrono>
#include <vector>
#include <map>
#include <sstream>


class AbstractKmerDb {
protected:

	uint32_t kmerLength;

	bool isInitialized;

	double fraction;

	double startFraction;

	std::vector<string> sampleNames;

	std::vector<size_t> sampleKmersCount;

	virtual void initialize(uint32_t kmerLength, double fraction) {
		this->kmerLength = kmerLength;
		this->fraction = fraction;
		this->isInitialized = true;
	}

public:

	AbstractKmerDb() : kmerLength(0), isInitialized(false), fraction(0), startFraction(0)  {}

	virtual ~AbstractKmerDb() {}

	uint32_t getKmerLength() const { return kmerLength; }

	double getFraction() const { return fraction; }

	double getStartFraction() const { return startFraction; }

	size_t getSamplesCount() const { return sampleNames.size(); }

	const std::vector<string>& getSampleNames() const { return sampleNames; }

	const std::vector<size_t>& getSampleKmersCount() const { return sampleKmersCount; }

	virtual size_t getKmersCount() const = 0;

	virtual size_t getPatternsCount() const = 0;

	virtual size_t getPatternBytes() const = 0;

	virtual size_t getHashtableBytes() const = 0;

	virtual size_t getHashtableEntrySize() const = 0;

	virtual void serialize(std::ofstream& file, bool rawHashtables) const = 0;

	virtual bool deserialize(std::ifstream& file, bool skipHashtables = false) = 0;

	virtual std::string printStats() const = 0;

	virtual std::string printDetailedTimes() const = 0;

	virtual std::string printProgress() const = 0;

	
	virtual sample_id_t addKmers(
		const std::string& sampleName,
		const kmer_t* kmers,
		size_t kmersCount,
		uint32_t kmerLength, 
		double fraction) {
		LOG_VERBOSE << "Adding sample " << sampleNames.size() + 1 << ": " << sampleName << " (" << kmersCount << " kmers)" << endl;
		
		if (!isInitialized) {
			initialize(kmerLength, fraction);
		}

		if (this->kmerLength != kmerLength) {
			throw std::runtime_error("Error in AbstractKmerDb::addKmers(): adding kmers of different length");
		}
		if (this->fraction != fraction) {
			throw std::runtime_error("Error in AbstractKmerDb::addKmers(): adding kmers of different minhash fraction");
		}

		sample_id_t newId = (sample_id_t) sampleNames.size();
		sampleNames.push_back(sampleName);
		sampleKmersCount.push_back(kmersCount);

		return newId;
	}
};
