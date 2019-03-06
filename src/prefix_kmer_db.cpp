#include "prefix_kmer_db.h"

#include "log.h"

#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif
#include <omp.h>

#include <numeric>


#define USE_PREFETCH
#define ALL_STATS

using namespace std;

const size_t PrefixKmerDb::ioBufferBytes = (2 << 29); //512MB buffer 


// *****************************************************************************************
//
PrefixKmerDb::PrefixKmerDb(int _num_threads) : 
	num_threads(_num_threads > 0 ? _num_threads : std::thread::hardware_concurrency())
{
	patterns.reserve(2 << 20);
	patterns.push_back(pattern_t());

	threadPatterns.resize(num_threads);
	for (auto& tp : threadPatterns) { 
		tp.reserve(2 << 20); 
	}
	
	workers.prefixHistogram.resize(num_threads);
	workers.hashtableSearch.resize(num_threads);
	workers.patternExtension.resize(num_threads);

	for (auto& t : workers.prefixHistogram) {
		t = std::thread(&PrefixKmerDb::histogramJob, this);
	}

	for (auto& t : workers.hashtableSearch) {
		t = std::thread(&PrefixKmerDb::hashtableSearchJob, this);
	}

	for (auto& t : workers.patternExtension) {
		t = std::thread(&PrefixKmerDb::patternJob, this);
	}
}

// *****************************************************************************************
//
PrefixKmerDb::~PrefixKmerDb() {
	queues.prefixHistogram.MarkCompleted();
	for (auto& t : workers.prefixHistogram) {
		t.join();
	}

	queues.hashtableSearch.MarkCompleted();
	for (auto& t : workers.hashtableSearch) {
		t.join();
	}

	queues.patternExtension.MarkCompleted();
	for (auto& t : workers.patternExtension) {
		t.join();
	}
}


// *****************************************************************************************
//
void PrefixKmerDb::initialize(uint32_t kmerLength, double fraction) {
	AbstractKmerDb::initialize(kmerLength, fraction);

	size_t prefixBits = (kmerLength - SUFFIX_LEN) * 2;
	size_t binsCount = 1 << prefixBits;
	this->prefixMask = ((1ULL << prefixBits) - 1) << SUFFIX_BITS;

	prefixHistogram.resize(binsCount);
	hashtables.resize(binsCount);
}

// *****************************************************************************************
//
void PrefixKmerDb::histogramJob() {
	while (!this->queues.prefixHistogram.IsCompleted()) {
		DictionarySearchTask task;

		if (this->queues.prefixHistogram.Pop(task)) {

			// generate histogram
			const std::vector<kmer_t>& kmers = *(task.kmers);

			size_t n_kmers = kmers.size();
			size_t block_size = n_kmers / task.num_blocks;
			size_t lo = task.block_id * block_size;
			size_t hi = (task.block_id == task.num_blocks - 1) ? n_kmers : lo + block_size;

			kmer_t boundaryPrefix = kmers[hi - 1] & this->prefixMask;
			uint32_t boundaryCount = 0;
			uint32_t begin = lo;

			kmer_t prevPrefix = kmers[lo] & this->prefixMask;
			
			for (size_t i = lo + 1; i < hi; ++i) {
				kmer_t prefix = kmers[i] & this->prefixMask;

				if (prefix != prevPrefix) {
					
					prefixHistogram[prevPrefix >> SUFFIX_BITS] = i - begin;
					begin = i;
					prevPrefix = prefix;
					
					if (prefix == boundaryPrefix) { // check if we are at the histogram bound 
						// all kmers starting from here have same prefix as the last one in the block
						boundaryCount = hi - i;
						break;
					}
				}
			}

			this->prefixHistogramMutex.lock();
			prefixHistogram[boundaryPrefix >> SUFFIX_BITS] += boundaryCount;
			this->prefixHistogramMutex.unlock();

			// resize hastables
			size_t tablesPerWorker = hashtables.size() / task.num_blocks;
			lo = task.block_id * tablesPerWorker;
			hi = (task.block_id == task.num_blocks - 1) ? hashtables.size() : lo + tablesPerWorker;

			this->mem.hashtable = 0;

			for (int i = lo; i < hi; ++i) {
				hashtables[i].reserve_for_additional(prefixHistogram[i]);
				this->mem.hashtable += hashtables[i].get_bytes();
			}
			
			this->semaphore.dec();
		}
	}
}

// *****************************************************************************************
//
void PrefixKmerDb::hashtableSearchJob() {
	while (!this->queues.hashtableSearch.IsCompleted()) {
		DictionarySearchTask task;

		if (this->queues.hashtableSearch.Pop(task)) {
			LOG_DEBUG << "Block " << task.block_id << " started" << endl;
			const std::vector<kmer_t>& kmers = *(task.kmers);

			size_t n_kmers = kmers.size();
			size_t block_size = n_kmers / task.num_blocks;
			size_t lo = task.block_id * block_size;
			size_t hi = (task.block_id == task.num_blocks - 1) ? n_kmers : lo + block_size;

			size_t existing_id = lo;
			size_t to_add_id = hi - 1;

			kmer_t u_kmer;
#ifdef USE_PREFETCH
			kmer_t prefetch_kmer;
			const size_t prefetch_dist = 48;
#endif

			for (size_t i = lo; i < hi; ++i)
			{
				u_kmer = kmers[i];
				

#ifdef USE_PREFETCH
				if (i + prefetch_dist < hi)
				{
					prefetch_kmer = kmers[i + prefetch_dist];
			//		kmers2patternIds.prefetch(prefetch_kmer);
				}
#endif
				// Check whether k-mer exists in a dictionary
				pattern_id_t* i_kmer;// = kmers2patternIds.find(u_kmer);
				pattern_id_t p_id;

				if (i_kmer == nullptr)
				{
					// do not add kmer to hashtable - just mark as to be added
					samplePatterns[to_add_id].first = u_kmer;
					--to_add_id;
				}
				else {
					p_id = *i_kmer;

					samplePatterns[existing_id].first = p_id;
					samplePatterns[existing_id].second = i_kmer;
					existing_id++;
				}
			}

			*(task.num_existing_kmers) = existing_id;
			//cout << "Thread " << tid << ", existing: " << existing_id - lo << ", to add: " << hi - existing_id << endl;

			LOG_DEBUG << "Block " << task.block_id << " finished" << endl;
			this->semaphore.dec();
		}
	}
}

// *****************************************************************************************
//
void PrefixKmerDb::patternJob() {
	while (!this->queues.patternExtension.IsCompleted()) {
		PatternExtensionTask task;

		if (this->queues.patternExtension.Pop(task)) {
			LOG_DEBUG << "Block " << task.block_id << " started" << endl;
			threadPatterns[task.block_id].clear();
			threadPatterns[task.block_id].reserve(task.ranges->back());

			size_t lo = (*task.ranges)[task.block_id];
			size_t hi = (*task.ranges)[task.block_id + 1];
			size_t mem = 0;

			for (size_t i = lo; i < hi;) {
				size_t j;
				auto p_id = samplePatterns[i].first;

				// count k-mers from current sample with considered template 
				for (j = i + 1; j < hi; ++j) {
					if (p_id != samplePatterns[j].first) {
						break;
					}
				}
				size_t pid_count = j - i;

				if (patterns[p_id].get_num_kmers() == pid_count && !patterns[p_id].get_is_parrent()) {
					// Extend pattern - all k-mers with considered template exist in the analyzed sample
					mem -= patterns[p_id].get_bytes();
					patterns[p_id].expand((sample_id_t)(task.sample_id));
					mem += patterns[p_id].get_bytes();
				}
				else
				{
					// Generate new template 
					pattern_id_t local_pid = task.new_pid->fetch_add(1);

					threadPatterns[task.block_id].emplace_back(local_pid, pattern_t(patterns[p_id], p_id, task.sample_id, (uint32_t)pid_count));
					mem += threadPatterns[task.block_id].back().second.get_bytes();

					if (p_id) {
						patterns[p_id].set_num_kmers(patterns[p_id].get_num_kmers() - pid_count);
					}

					for (size_t k = i; k < j; ++k) {
						*(samplePatterns[k].second) = local_pid;
					}

				}

				i = j;
			}

			(*task.threadBytes)[task.block_id] = mem;

			LOG_DEBUG << "Block " << task.block_id << " finished" << endl;
			this->semaphore.dec();
		}
	}
}


// *****************************************************************************************
//
sample_id_t PrefixKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers, kmerLength, fraction);

	// get prefix histogram (parallel)
	LOG_DEBUG << "Restructurizing hashtable (parallel)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	// prepare tasks
	for (int tid = 0; tid < num_threads; ++tid) {
		semaphore.inc();
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.prefixHistogram.Push(DictionarySearchTask{ tid, num_threads, &kmers, nullptr });
	}

	semaphore.waitForZero();
	times.hashtableResize += std::chrono::high_resolution_clock::now() - start;

#ifdef _DEBUG
	uint32_t histoSum = std::accumulate(prefixHistogram.begin(), prefixHistogram.end(), 0);
	if (histoSum != kmers.size()) {
		throw std::runtime_error("PrefixKmerDb::addKmers() - invalid histogram sum");
	}
#endif


	return 0;
}



// *****************************************************************************************
//
void PrefixKmerDb::serialize(std::ofstream& file) const {

	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<kmer_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());

	// store number of samples
	size_t temp = getSamplesCount();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

	// store sample info 
	for (int i = 0; i < sampleNames.size(); ++i) {
		temp = sampleKmersCount[i]; // store kmer count
		file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

		const string& s = sampleNames[i]; // store name
		temp = s.size();
		file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));
		file.write(s.data(), temp);
	}

	// store number of hashmaps
	temp = hashtables.size();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

	// store all hashmaps
	for (const auto& ht : hashtables) {

		// store ht size
		temp = ht.get_size();
		file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

		// write ht elements in portions
		size_t bufpos = 0;
		for (auto it = ht.cbegin(); it < ht.cend(); ++it) {
			if (ht.is_free(*it)) {
				continue;
			}

			hashtableBuffer[bufpos++] = *it;
			if (bufpos == numHastableElements) {
				file.write(reinterpret_cast<const char*>(&bufpos), sizeof(size_t));
				file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
				bufpos = 0;
			}
		}
		// write remaining ht elements
		file.write(reinterpret_cast<const char*>(&bufpos), sizeof(size_t));
		file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
	}

	// write patterns in portions
	temp = patterns.size();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

	char * currentPtr = buffer;
	for (int pid = 0; pid < patterns.size(); ++pid) {
		if (currentPtr + patterns[pid].get_bytes() > buffer + ioBufferBytes) {
			size_t blockSize = currentPtr - buffer;
			file.write(reinterpret_cast<const char*>(&blockSize), sizeof(size_t)); // write size of block to facilitate deserialization
			file.write(buffer, blockSize);
			currentPtr = buffer;
		}

		currentPtr = patterns[pid].pack(currentPtr);
		// this should never happen
		if (currentPtr > buffer + ioBufferBytes) {
			throw std::runtime_error("Buffer overflow when saving patterns!");
		}
	}

	// write remaining patterns
	size_t blockSize = currentPtr - buffer;
	file.write(reinterpret_cast<const char*>(&blockSize), sizeof(size_t)); // write size of block to facilitate deserialization
	file.write(buffer, blockSize);

	// save kmer length and fraction
	file.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength)); // write size of block to facilitate deserialization

	file.write(reinterpret_cast<const char*>(&fraction), sizeof(fraction)); // write size of block to facilitate deserialization
}

// *****************************************************************************************
//
bool PrefixKmerDb::deserialize(std::ifstream& file) {

	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<kmer_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());

	// load sample info
	size_t temp;
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	sampleNames.resize(temp);
	sampleKmersCount.resize(temp);

	for (int i = 0; i < sampleNames.size(); ++i) {
		file.read(reinterpret_cast<char*>(&temp), sizeof(temp));  // load kmer count
		sampleKmersCount[i] = temp;

		string& s = sampleNames[i];
		file.read(reinterpret_cast<char*>(&temp), sizeof(temp));  // load sample name
		file.read(buffer, temp);
		s.assign(buffer, temp);
	}

	if (!file) {
		return false;
	}

	// load number of hashmaps
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	hashtables.resize(temp);

	// load all hashtables
	for (auto& ht : hashtables) {

		// loal ht size
		file.read(reinterpret_cast<char*>(&temp), sizeof(temp));

		// load ht elements
		ht.clear();
		ht.reserve_for_additional(temp);

		size_t readCount = 0;
		while (readCount < temp) {
			size_t portion = 0;
			file.read(reinterpret_cast<char*>(&portion), sizeof(size_t));
			file.read(buffer, portion * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));

			for (size_t j = 0; j < portion; ++j) {
				ht.insert(hashtableBuffer[j].key, hashtableBuffer[j].val);
			}
			readCount += portion;
		}
	}

	if (!file) {
		return false;
	}

	// load patterns
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	patterns.clear();
	patterns.resize(temp);

	size_t pid = 0;
	while (pid < patterns.size()) {
		size_t blockSize;
		file.read(reinterpret_cast<char*>(&blockSize), sizeof(size_t));
		file.read(buffer, blockSize);

		char * currentPtr = buffer;
		while (currentPtr < buffer + blockSize) {
			currentPtr = patterns[pid].unpack(currentPtr);
			++pid;
		}
	}

	if (!file) {
		return false;
	}

	// load kmer length and fraction
	file.read(reinterpret_cast<char*>(&kmerLength), sizeof(kmerLength));

	file.read(reinterpret_cast<char*>(&fraction), sizeof(fraction));

	return true;
}

