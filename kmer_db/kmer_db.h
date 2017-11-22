#pragma once
#include "pattern.h"
#include "hashmap_lp.h"
#include "hashmap_dh.h"
#include "array.h"
#include "queue.h"
#include "aligned_vector.h"
#include "parallel_sorter.h"

#include <map>
#include <fstream>
#include <chrono>
#include <vector>
#include <map>



typedef uint64_t kmer_t;


class AbstractKmerDb {
protected:
	std::vector<string> sampleNames;

public:

	AbstractKmerDb() {}

	const size_t getSamplesCount() const { return sampleNames.size(); }

	const std::vector<string>& getSampleNames() const { return sampleNames; }

	virtual const size_t getKmersCount() const = 0;

	virtual const size_t getPatternsCount() const = 0;

	virtual const size_t getPatternBytes() const = 0;

	virtual const size_t getHashtableBytes() const = 0;

	virtual std::vector<kmer_t> getKmers() const = 0;

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers) {
		sample_id_t newId = sampleNames.size();
		sampleNames.push_back(sampleName);
		return newId;

	}

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const = 0;

	//virtual void calculateSimilarity(Array<uint32_t>& matrix) const = 0;
};


class NaiveKmerDb : public AbstractKmerDb {
public:
	NaiveKmerDb() : kmers2patternIds((unsigned long long) - 1) {
		patternBytes = 0;
	}

	virtual const size_t getKmersCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternsCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternBytes() const { return patternBytes; }

	virtual const size_t getHashtableBytes() const { return kmers2patternIds.get_bytes(); }

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

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarity(Array<uint32_t>& matrix) const;

	std::map<std::vector<sample_id_t>, size_t> getPatternsStatistics() const;

protected:
	size_t patternBytes;		// iloœæ pamiêci zajmowana przez wszystkie wzorce

	hash_map_lp<kmer_t, pattern_id_t> kmers2patternIds;

	std::vector<std::vector<sample_id_t>> patterns;


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

	FastKmerDb();

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

	virtual sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix) const;
	virtual void calculateSimilarityDirect(LowerTriangularMatrix<uint32_t>& matrix) const;
	virtual void calculateSimilarityBuffered(LowerTriangularMatrix<uint32_t>& matrix) const;

	virtual void calculateSimilarity(const FastKmerDb& sampleDb, std::vector<uint32_t>& vector) const;
	virtual void calculateSimilarityOld(const FastKmerDb& sampleDb, std::vector<uint32_t>& vector) const;


	virtual void serialize(std::ofstream& file) const;

	virtual bool deserialize(std::ifstream& file);

	virtual void savePatterns(std::ofstream& file) const;

protected:
	static const size_t ioBufferBytes;
	
	size_t patternBytes;		// iloœæ pamiêci zajmowana przez wBuszystkie wzorce
								// K-mer database structures
	hash_map_lp<kmer_t, pattern_id_t> kmers2patternIds;

	// liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
	std::vector<pattern_t> patterns;

	std::vector<std::vector<std::pair<pattern_id_t, pattern_t>>> threadPatterns;

	// first element - pattern id, second element 
//	std::vector<std::pair<pattern_id_t, pattern_id_t*>> samplePatterns;
	aligned_vector<std::pair<pattern_id_t, pattern_id_t*>> samplePatterns;
#ifdef USE_RADULS
	aligned_vector<std::pair<pattern_id_t, pattern_id_t*>> tmp_samplePatterns;
#endif

	CRegisteringQueue<DictionarySearchTask> dictionarySearchQueue;
	
	CRegisteringQueue<PatternExtensionTask> patternExtensionQueue;

	std::vector<std::thread> dictionarySearchWorkers;

	std::vector<std::thread> patternExtensionWorkers;

	std::vector<int> kmers_to_add_to_HT;

	Semaphore semaphore;

	int num_threads;
};