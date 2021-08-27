#pragma once
#include "kmer_db.h"

#include <atomic>
#include <numeric>

// auxiliary templates for saving/loading single values in a binary file
template <class T>
void save(std::ostream &stream, const T& val) {
	stream.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <class T>
void load(std::istream &stream, T& val) {
	stream.read(reinterpret_cast<char*>(&val), sizeof(T));
}

struct HashtableTask {
	uint32_t block_id;
	uint32_t lo;
	uint32_t hi;
	const kmer_t* kmers;
	size_t n_kmers;
};

struct PatternTask {
	int block_id;
	uint32_t lo;
	uint32_t hi;
	sample_id_t sample_id;
	std::atomic<pattern_id_t>* new_pid;
};


class PrefixKmerDb : public AbstractKmerDb {
public:
	size_t getKmersCount() const override { return kmersCount; }

	size_t getPatternsCount() const override { return patterns.size(); }

	size_t getPatternBytes() const override { return stats.patternBytes; }

	size_t getHashtableBytes() const override { return stats.hashtableBytes; };

	size_t getWorkersPatternBytes() const {
		size_t total = 0;
		for (const auto& t : threadPatterns) {
			total += t.capacity() * sizeof(std::pair<pattern_id_t, pattern_t>);
		}
		return total;

	}

	size_t getHashtableEntrySize() const override { return sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t); };
	
	const std::vector<hash_map_lp<suffix_t, pattern_id_t>>& getHashtables() const { return hashtables; }

	const std::vector<pattern_t>& getPatterns() const { return patterns; }

	std::vector<pattern_t>& getPatterns() { return patterns; }

	PrefixKmerDb(int _num_threads);

	~PrefixKmerDb();

	sample_id_t addKmers(const std::string& sampleName,
		const kmer_t* kmers,
		size_t kmersCount,
		uint32_t kmerLength,
		double fraction) override;

	void serialize(std::ofstream& file, bool rawHashtables) const override;

	bool deserialize(std::ifstream& file, bool skipHashtables = false) override;

	std::string printProgress() const override {
		std::ostringstream oss;

		size_t otherBytes =
			prefixHistogram.capacity() * sizeof(uint32_t) +
			hashtables.capacity() * sizeof(hash_map_lp<suffix_t, pattern_id_t>) +
			samplePatterns.get_bytes();
	
		oss << "HT entries: " << Log::formatLargeNumber(getKmersCount())
			<< " (" << ((getKmersCount() * getHashtableEntrySize()) >> 20) << " MB, " << (getHashtableBytes() >> 20) << " MB res),"
			<< "\t patterns: " << Log::formatLargeNumber(getPatternsCount())
			<< " (" << Log::formatLargeNumber(getPatternBytes()) << " B), worker_patterns: " << (getWorkersPatternBytes() >> 20) << " MB res, "
			<< " other: " << (otherBytes >> 20) << " MB";
		return oss.str();
	}

	std::string printStats() const override {
		std::ostringstream oss;
		oss << "Number of samples: " << Log::formatLargeNumber(getSamplesCount()) << endl
			<< "Number of patterns: " << Log::formatLargeNumber(getPatternsCount()) << " (" << Log::formatLargeNumber(stats.patternBytes) << " B)" << endl
			<< "Number of k-mers: " << Log::formatLargeNumber(getKmersCount()) << endl
			<< "K-mer length: " << Log::formatLargeNumber(getKmerLength()) << endl
			<< "Minhash fraction: " << getFraction() << endl
			<< "Workers count: " << num_threads << endl;
		return oss.str();
	}

	std::string printDetailedTimes() const override {
		std::ostringstream oss;
		oss << "\tHashatable processing (parallel): " << times.hashtableProcess.count() <<  endl 
		//	<< "\timbalance: " << stats.hashtableJobsImbalance / getSamplesCount() <<  endl
			<< "\t\tResize: " << (times.hashtableResize_worker.count() / num_threads) << endl
			<< "\t\tFind'n'add: " << (times.hashtableFind_worker.count() / num_threads) << endl
			<< "\tSort time (parallel): " << times.sort.count() << endl
			<< "\tPattern extension time (parallel): " << times.extension.count() << endl;
		return oss.str();
	}

	void savePatterns(std::ofstream& file) const;

protected:

	static const size_t IO_BUFFER_BYTES = (2 << 29); //512MB buffer 
	static const int PREFETCH_DIST = 48;

	static const uint64_t SERIALIZATION_RAW_HASHTABLES = 0x01;

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
		RegisteringQueue<HashtableTask> hashtableAddition{ 1 };

		RegisteringQueue<PatternTask> patternExtension{ 1 };
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
		
		std::chrono::duration<double> hashtableResize_worker{ 0 };
		std::chrono::duration<double> hashtableFind_worker{ 0 };
		std::chrono::duration<double> hashtableAdd_worker{ 0 };


		std::chrono::duration<double> sort { 0 };
		std::chrono::duration<double> extension{ 0 };
	} times;

	// structure for storing bytes 
	struct {
		std::atomic<size_t> hashtableBytes { 0 };
		std::atomic<size_t> patternBytes { 0 };

		double hashtableJobsImbalance{ 0 };

	} stats;


	void initialize(uint32_t kmerLength, double fraction) override;

	//void histogramJob();
	//void hashtableSearchJob();
	//void hashtableAdditionJob();

	void hashtableJob();
	void patternJob();

};
