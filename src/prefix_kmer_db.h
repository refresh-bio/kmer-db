#pragma once
#include "kmer_db.h"

#define SUFFIX_LEN 16
#define SUFFIX_BITS (SUFFIX_LEN * 2)




class PrefixKmerDb : public AbstractKmerDb {
public:

	const size_t getKmersCount() const override { return kmersCount; }

	const size_t getPatternsCount() const override { return patterns.size(); }

	const size_t getPatternBytes() const override { return mem.pattern; }

	const size_t getHashtableBytes() const override { return mem.hashtable; };
	
	const std::vector<hash_map_lp<kmer_t, pattern_id_t>>& getHashtables() const { return hashtables; }

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
		oss << "\tHashatable resizing (parallel): " << times.hashtableResize.count() << endl
			<< "\tHashtable searching (parallel): " << times.hashtableFind.count() << endl
			<< "\tHashatable insertion (parallel): " << times.hashtableAdd.count() << endl
			<< "\tSort time (parallel): " << times.sort.count() << endl
			<< "\tPattern extension time (parallel): " << times.extension.count() << endl;
		return oss.str();
	}

protected:

	static const size_t ioBufferBytes;

	int num_threads;
	
	size_t kmersCount;

	uint64_t prefixMask;

	std::vector<uint32_t> prefixHistogram;

	std::vector<hash_map_lp<kmer_t, pattern_id_t>> hashtables;

	std::vector<pattern_t> patterns;

	std::vector<std::vector<std::pair<pattern_id_t, pattern_t>>> threadPatterns;

	// first element - pattern id, second element - ht table entry
	aligned_vector<std::pair<pattern_id_t, pattern_id_t*>> samplePatterns;

	// struct for storing queues
	struct {
		CRegisteringQueue<DictionarySearchTask> prefixHistogram{ 1 };

		CRegisteringQueue<DictionarySearchTask> hashtableSearch{ 1 };

		CRegisteringQueue<PatternExtensionTask> patternExtension{ 1 };
	} queues;

	// struct for storing workers
	struct {
		std::vector<std::thread> prefixHistogram;

		std::vector<std::thread> hashtableSearch;

		std::vector<std::thread> patternExtension;
	} workers;

	std::vector<int> kmers_to_add_to_HT;

	std::mutex prefixHistogramMutex;

	Semaphore semaphore;

	// structure for storing all the times
	struct {
		std::chrono::duration<double> hashtableResize { 0 };
		std::chrono::duration<double> hashtableFind { 0 };
		std::chrono::duration<double> hashtableAdd { 0 };
		std::chrono::duration<double> sort { 0 };
		std::chrono::duration<double> extension{ 0 };
	} times;

	// structure for storing bytes 
	struct {
		size_t hashtable { 0 };
		size_t pattern { 0 };
	} mem;


	void initialize(uint32_t kmerLength, double fraction) override;

	void histogramJob();
	void hashtableSearchJob();
	void hashtableAdditionJob();
	void patternJob();
};
