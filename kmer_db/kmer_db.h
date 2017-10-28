#pragma once
#include "pattern.h"
#include "hashmap_lp.h"
#include "hashmap_dh.h"
#include "Array.h"

#include <map>

typedef uint16_t sample_id_t;
typedef uint64_t kmer_t;
typedef uint64_t pattern_id_t;


class AbstractKmerDb {
protected:
	size_t numSamples;

public:

	AbstractKmerDb() : numSamples(0) {}

	const size_t getNumSamples() const { return numSamples; }

	virtual const size_t getKmersCount() const = 0;

	virtual const size_t getPatternsCount() const = 0;

	virtual const size_t getPatternMem() const = 0;

	virtual const size_t getHashtableMem() const = 0;

	virtual bool extractKmers(const std::string &filename, std::vector<kmer_t>& kmers);

	virtual void addKmers(sample_id_t sampleId, const std::vector<kmer_t>& kmers) { ++numSamples; }

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const = 0;

	virtual void calculateSimilarityMatrix(Array<uint32_t>& matrix) const = 0;
};


class NaiveKmerDb : public AbstractKmerDb {
public:
	size_t mem_pattern_desc = 0;		// iloœæ pamiêci zajmowana przez wszystkie wzorce
	
	hash_map<kmer_t, pattern_id_t> kmers2patternIds;

	std::vector<std::vector<sample_id_t>> patterns;

public:
	NaiveKmerDb() : kmers2patternIds((unsigned long long) - 1, (unsigned long long) - 2) {
		mem_pattern_desc = 0;
	}

	virtual const size_t getKmersCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternsCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternMem() const { return mem_pattern_desc; }

	virtual const size_t getHashtableMem() const { return kmers2patternIds.getMem(); }

	virtual void addKmers(sample_id_t sampleId, const std::vector<kmer_t>& kmers);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarityMatrix(Array<uint32_t>& matrix) const;

	std::map<std::vector<sample_id_t>, size_t> getPatternsStatistics() const;
};


class FastKmerDb : public AbstractKmerDb {
protected:
	size_t mem_pattern_desc;		// iloœæ pamiêci zajmowana przez wszystkie wzorce
	
	// K-mer database structures
	hash_map_dh<kmer_t, pattern_id_t> kmers2patternIds;

	// liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
	std::vector<pattern_t<sample_id_t>> patterns;		

	std::vector<std::vector<std::pair<pattern_id_t, pattern_t<sample_id_t>>>> threadPatterns;

	// first element - pattern id, second element 
	std::vector<std::pair<pattern_id_t, pattern_id_t*>> samplePatterns;

	
public:

	FastKmerDb();

	virtual const size_t getKmersCount() const { return kmers2patternIds.get_size(); }

	virtual const size_t getPatternsCount() const { return patterns.size(); }

	virtual const size_t getPatternMem() const { return mem_pattern_desc; }

	virtual const size_t getHashtableMem() const { return kmers2patternIds.getMem(); };

	virtual void addKmers(sample_id_t sampleId, const std::vector<kmer_t>& kmers);

	virtual void mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const;

	virtual void calculateSimilarityMatrix(Array<uint32_t>& matrix) const;

	std::map<std::vector<sample_id_t>, size_t> getPatternsStatistics() const;

};