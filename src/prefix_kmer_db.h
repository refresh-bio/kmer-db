#pragma once
#include "kmer_db.h"

#include <atomic>

struct HashtableAdditionTask {
	int block_id;
	const std::vector<kmer_t>* kmers;
	const std::vector<size_t>* ranges;
};


class PrefixKmerDb : public AbstractKmerDb {
public:

	size_t getKmersCount() const override { return kmersCount; }

	size_t getPatternsCount() const override { return patterns.size(); }

	size_t getPatternBytes() const override { return mem.pattern; }

	size_t getHashtableBytes() const override { return mem.hashtable; };

	size_t getHashtableEntrySize() const override { return sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t); };
	
	const std::vector<hash_map_lp<suffix_t, pattern_id_t>>& getHashtables() const { return hashtables; }

	const std::vector<pattern_t>& getPatterns() const { return patterns; }

	std::vector<pattern_t>& getPatterns() { return patterns; }

	PrefixKmerDb(int _num_threads);

	~PrefixKmerDb();

	sample_id_t addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) override;

	void serialize(std::ofstream& file) const override;

	bool deserialize(std::ifstream& file) override;

	std::string printStats() const override {
		std::ostringstream oss;
		oss << "Number of samples: " << Log::formatLargeNumber(getSamplesCount()) << endl
			<< "Number of patterns: " << Log::formatLargeNumber(getPatternsCount()) << endl
			<< "Number of k-mers: " << Log::formatLargeNumber(getKmersCount()) << endl
			<< "K-mer length: " << Log::formatLargeNumber(getKmerLength()) << endl
			<< "Minhash fraction: " << getFraction() << endl;
		return oss.str();
	}

	std::string printDetailedTimes() const override {
		std::ostringstream oss;
		oss << "\tHashatable processing (parallel): " << times.hashtableProcess.count() << endl
			<< "\tSort time (parallel): " << times.sort.count() << endl
			<< "\tPattern extension time (parallel): " << times.extension.count() << endl;
		return oss.str();
	}

protected:

	static const size_t ioBufferBytes = (2 << 29); //512MB buffer 
	static const int prefetch_dist = 48;

	int num_threads;
	
	std::atomic<size_t> kmersCount{ 0 };

	std::vector<uint32_t> prefixHistogram;

	std::vector<hash_map_lp<suffix_t, pattern_id_t>> hashtables;

	std::vector<pattern_t> patterns;

	std::vector<std::vector<std::pair<pattern_id_t, pattern_t>>> threadPatterns;

	// first element - pattern id, second element - ht table entry
	aligned_vector<std::pair<kmer_or_pattern_t, pattern_id_t*>> samplePatterns;

	// struct for storing queues
	struct {
		CRegisteringQueue<HashtableAdditionTask> hashtableAddition{ 1 };

		CRegisteringQueue<PatternExtensionTask> patternExtension{ 1 };
	} queues;

	// struct for storing workers
	struct {
		std::vector<std::thread> hashtableAddition;

		std::vector<std::thread> patternExtension;
	} workers;

	Semaphore semaphore;

	// structure for storing all the times
	struct {
		std::chrono::duration<double> hashtableProcess { 0 };
		std::chrono::duration<double> sort { 0 };
		std::chrono::duration<double> extension{ 0 };
	} times;

	// structure for storing bytes 
	struct {
		std::atomic<size_t> hashtable { 0 };
		std::atomic<size_t> pattern { 0 };
	} mem;


	void initialize(uint32_t kmerLength, double fraction) override;

	void histogramJob();
	void hashtableSearchJob();
	void hashtableAdditionJob();

	void hashtableJob();
	void patternJob();
};
