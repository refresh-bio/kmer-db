#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "pattern.h"
#include "hashmap_lp.h"
//#include "hashset_lp.h"
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



#define KMER_MSB (1ULL << 63)

typedef uint64_t kmer_t;


class AbstractKmerDb {
protected:

	uint32_t kmerLength;

	double fraction;

	std::vector<string> sampleNames;

	std::vector<size_t> sampleKmersCount;

public:

	AbstractKmerDb() : kmerLength(0) {}

	const uint32_t getKmerLength() const { return kmerLength; }

	const double getFraction() const { return fraction; }

	const size_t getSamplesCount() const { return sampleNames.size(); }

	const std::vector<string>& getSampleNames() const { return sampleNames; }

	const std::vector<size_t>& getSampleKmersCount() const { return sampleKmersCount; }

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) {
		LOG_VERBOSE << "Adding sample: " << sampleName << " (" << kmers.size() << " kmers)" << endl;
		sample_id_t newId = (sample_id_t) sampleNames.size();
		sampleNames.push_back(sampleName);
		sampleKmersCount.push_back(kmers.size());

		if (this->kmerLength == 0) {
			this->kmerLength = kmerLength;
			this->fraction = fraction;
		}
		else {
			if (this->kmerLength != kmerLength) {
				throw std::runtime_error("ERROR in FastKmerDb::addKmers(): adding kmers of different length");
			}
			if (this->fraction != fraction) {
				throw std::runtime_error("ERROR in FastKmerDb::addKmers(): adding kmers of different minhash fraction");
			}
		}
		
		return newId;
	}

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const = 0;
};


struct DictionarySearchTask {
	int block_id;
	int num_blocks;
	const std::vector<kmer_t>* kmers;
	size_t* num_existing_kmers;
};

struct PatternExtensionTask {
	int block_id;
	size_t sample_id;
	std::vector<size_t>* ranges;
	std::atomic<size_t>* new_pid;
	std::vector<size_t>* threadBytes;

};

class FastKmerDb : public AbstractKmerDb {
public:
	FastKmerDb(int _num_threads, size_t cacheBufferMb);

	~FastKmerDb();

	const size_t getKmersCount() const { return kmers2patternIds.get_size(); }

	const size_t getPatternsCount() const { return patterns.size(); }

	const size_t getPatternBytes() const { return patternBytes; }

	const size_t getHashtableBytes() const { return kmers2patternIds.get_bytes(); };

	const size_t getRepeatedKmersCount() const { return 0; }

//	const hash_set_lp<kmer_t>& getRepeatedKmers() const { return repeatedKmers; }

	const hash_map_lp<kmer_t, pattern_id_t>& getKmers2patternIds() const { return kmers2patternIds; }

	const std::vector<pattern_t>& getPatterns() const { return patterns; }

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix);// const;

	virtual void calculateSimilarity(const std::vector<kmer_t>& kmers, std::vector<uint32_t>& vector) const;
	
	virtual void serialize(std::ofstream& file) const;

	virtual bool deserialize(std::ifstream& file);

	virtual void savePatterns(std::ofstream& file) const;

	std::vector<kmer_t> extractKmers() const {
		std::vector<kmer_t> kmers(kmers2patternIds.get_size());
		size_t i = 0;
		for (auto it = kmers2patternIds.cbegin(); it < kmers2patternIds.cend(); ++it) {
			if (kmers2patternIds.is_free(*it)) {
				continue;
			}
			kmers[i++] = it->key;
		}
		return kmers;
	}

	std::string printStats() const {
		std::ostringstream oss;
		oss << "Number of samples: " << Log::formatLargeNumber(getSamplesCount()) << endl
			<< "Number of patterns: " << Log::formatLargeNumber(getPatternsCount()) << endl
			<< "Number of k-mers: " << Log::formatLargeNumber(getKmersCount()) << endl
			<< "Number of sample-repeated k-mers: " << Log::formatLargeNumber(getRepeatedKmersCount()) << endl
			<< "K-mer length: " << Log::formatLargeNumber(getKmerLength()) << endl
			<< "Minhash fraction: " << getFraction() << endl;
		return oss.str();
	}

	std::string printDetailedTimes() const {
		std::ostringstream oss;
		oss << "\tHashatable resizing (serial): " << hashtableResizeTime.count() << endl
			<< "\tHashtable searching (parallel): " << hashtableFindTime.count() << endl
			<< "\tHashatable insertion (serial): " << hashtableAddTime.count() << endl
			<< "\tSort time (parallel): " << sortTime.count() << endl
			<< "\tPattern extension time (serial): " << extensionTime.count() << endl;
		return oss.str();
	}

protected:
	
	std::chrono::duration<double> hashtableResizeTime;
	std::chrono::duration<double> hashtableFindTime;
	std::chrono::duration<double> hashtableAddTime;
	std::chrono::duration<double> sortTime;
	std::chrono::duration<double> extensionTime;

	bool avx2_present;

	static const size_t ioBufferBytes;
	
	// memory needed for all templates
	size_t patternBytes;
								
	hash_map_lp<kmer_t, pattern_id_t> kmers2patternIds;

	//hash_set_lp<kmer_t> repeatedKmers;

	//std::vector<std::vector<kmer_t>> threadRepeatedKmers;

	// liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
	std::vector<pattern_t> patterns;

	std::vector<std::vector<std::pair<pattern_id_t, pattern_t>>> threadPatterns;

	// first element - pattern id, second element 
	aligned_vector<std::pair<pattern_id_t, pattern_id_t*>> samplePatterns;

	CRegisteringQueue<DictionarySearchTask> dictionarySearchQueue;
	
	CRegisteringQueue<PatternExtensionTask> patternExtensionQueue;

	std::vector<std::thread> dictionarySearchWorkers;

	std::vector<std::thread> patternExtensionWorkers;

	std::vector<int> kmers_to_add_to_HT;

	Semaphore semaphore;

	int num_threads;

	size_t cacheBufferMb;
};