#pragma once
#include "pattern.h"
#include "hashmap_lp.h"

typedef uint16_t sample_id_t;

typedef struct
{
	uint32_t num_occ;
	subpattern_t<sample_id_t> *last_subpattern;
} pattern_t;


class AbstractKmerDb {
public:

	virtual bool extractKmers(const std::string &filename, std::vector<uint64_t>& kmers);

	virtual void addKmers(sample_id_t sampleId, const std::vector<uint64_t>& kmers) = 0;

	virtual void getKmerSamples(uint64_t kmer, std::vector<sample_id_t>& samples) = 0;

	virtual const size_t getKmersCount() const = 0;
};


class NaiveKmerDb : public AbstractKmerDb {
public:
	std::unordered_map<uint64_t, std::vector<sample_id_t>> kmers2samples;

public:
	virtual void addKmers(sample_id_t sampleId, const std::vector<uint64_t>& kmers);

	virtual void getKmerSamples(uint64_t kmer, std::vector<sample_id_t>& samples);

	virtual const size_t getKmersCount() const { return kmers2samples.size(); }


};


class FastKmerDb : public AbstractKmerDb {
protected:
	// K-mer database structures
	hash_map<uint64_t, uint32_t> hm_kmer_dict;

	// liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
	std::vector<pattern_t> v_kmer_patterns;					

	std::vector<std::pair<uint64_t, uint32_t*>> v_current_file_pids;

	
public:

	FastKmerDb();

	virtual void addKmers(sample_id_t sampleId, const std::vector<uint64_t>& kmers);

	virtual void getKmerSamples(uint64_t kmer, std::vector<sample_id_t>& samples);

	virtual const size_t getKmersCount() const { return hm_kmer_dict.get_size(); }
	virtual const size_t getPatternsCount() const { return v_kmer_patterns.size(); }


};