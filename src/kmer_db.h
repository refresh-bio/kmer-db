#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

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

	virtual const size_t getKmersCount() const = 0;

	virtual const size_t getPatternsCount() const = 0;

	virtual const size_t getPatternBytes() const = 0;

	virtual const size_t getHashtableBytes() const = 0;

	virtual std::vector<kmer_t> getKmers() const = 0;

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) {
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
	std::vector<size_t>* num_existing_kmers;

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

	std::chrono::duration<double> hashtableResizeTime;
	std::chrono::duration<double> hashtableFindTime;
	std::chrono::duration<double> hashtableAddTime;
	std::chrono::duration<double> sortTime;
	std::chrono::duration<double> extensionTime;

	bool avx2_present;

	FastKmerDb(int _num_threads, size_t cacheBufferMb);

	~FastKmerDb();

	virtual const size_t getKmersCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternsCount() const { return patterns.size(); }

	virtual const size_t getPatternBytes() const { return patternBytes; }

	virtual const size_t getHashtableBytes() const { return kmers2patternIds.get_bytes(); };

	virtual std::vector<kmer_t> getKmers() const {
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

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix);// const;

	virtual void calculateSimilarity(const FastKmerDb& sampleDb, std::vector<uint32_t>& vector) const;
	
	virtual void serialize(std::ofstream& file) const;

	virtual bool deserialize(std::ifstream& file);

	virtual void savePatterns(std::ofstream& file) const;

protected:
	static const size_t ioBufferBytes;
	
	// memory needed for all templates
	size_t patternBytes;
								
	hash_map_lp<kmer_t, pattern_id_t> kmers2patternIds;

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