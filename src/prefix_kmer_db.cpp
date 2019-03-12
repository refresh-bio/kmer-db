#include "prefix_kmer_db.h"

#include "log.h"

#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif
#include <omp.h>

#include <numeric>
#include <cassert>


#define USE_PREFETCH
#define ALL_STATS

using namespace std;

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
	
//	workers.prefixHistogram.resize(num_threads);
//	workers.hashtableSearch.resize(num_threads);
	workers.hashtableAddition.resize(num_threads);
	workers.patternExtension.resize(num_threads);
/*
	for (auto& t : workers.prefixHistogram) {
		t = std::thread(&PrefixKmerDb::histogramJob, this);
	}

	for (auto& t : workers.hashtableSearch) {
		t = std::thread(&PrefixKmerDb::hashtableSearchJob, this);
	}

	for (auto& t : workers.hashtableAddition) {
		t = std::thread(&PrefixKmerDb::hashtableAdditionJob, this);
	}
	*/

	for (auto& t : workers.hashtableAddition) {
		t = std::thread(&PrefixKmerDb::hashtableJob, this);
	}

	for (auto& t : workers.patternExtension) {
		t = std::thread(&PrefixKmerDb::patternJob, this);
	}
}

// *****************************************************************************************
//
PrefixKmerDb::~PrefixKmerDb() {
/*
	queues.prefixHistogram.MarkCompleted();
	for (auto& t : workers.prefixHistogram) {
		t.join();
	}

	queues.hashtableSearch.MarkCompleted();
	for (auto& t : workers.hashtableSearch) {
		t.join();
	}
*/
	queues.hashtableAddition.MarkCompleted();
	for (auto& t : workers.hashtableAddition) {
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

	prefixHistogram.resize(binsCount);
	hashtables.resize(binsCount);
}

// *****************************************************************************************
//
void PrefixKmerDb::histogramJob() {
	while (!this->queues.prefixHistogram.IsCompleted()) {
		DictionarySearchTask task;

		if (this->queues.prefixHistogram.Pop(task)) {

			// determine ranges
			const std::vector<kmer_t>& kmers = *(task.kmers);

			size_t block_size = kmers.size() / task.num_blocks;
			size_t lo = task.block_id * block_size;
			size_t hi = (task.block_id == task.num_blocks - 1) ? kmers.size() : lo + block_size;

			// generate histogram
			kmer_t boundaryPrefix = GET_PREFIX_SHIFTED(kmers[hi - 1]);
			uint32_t boundaryCount = 0;
			uint32_t begin = lo;

			kmer_t prevPrefix = GET_PREFIX_SHIFTED(kmers[lo]);

			// only one prefix in a block
			if (prevPrefix == boundaryPrefix) {
				boundaryCount = hi - lo;
			}
			else {
				for (size_t i = lo + 1; i < hi; ++i) {
					kmer_t prefix = GET_PREFIX_SHIFTED(kmers[i]);

					if (prefix != prevPrefix) {

						prefixHistogram[prevPrefix] = i - begin;
						begin = i;
						prevPrefix = prefix;

						if (prefix == boundaryPrefix) { // check if we are at the histogram bound 
							// all kmers starting from here have same prefix as the last one in the block
							boundaryCount = hi - i;
							break;
						}
					}
				}
			}

			internalSempahores[0].dec_notify_all();
			internalSempahores[0].waitForZero();

			this->internalMutex.lock();
			prefixHistogram[boundaryPrefix] += boundaryCount;
			this->internalMutex.unlock();

			internalSempahores[1].dec_notify_all();
			internalSempahores[1].waitForZero();

			// resize hastables
			size_t tablesPerWorker = hashtables.size() / task.num_blocks;
			lo = task.block_id * tablesPerWorker;
			hi = (task.block_id == task.num_blocks - 1) ? hashtables.size() : lo + tablesPerWorker;

			size_t htBytes = 0;

			for (int i = lo; i < hi; ++i) {
				hashtables[i].reserve_for_additional(prefixHistogram[i]);
				htBytes += hashtables[i].get_bytes();
			}

			// update memory statistics (atomic - no sync needed)
			mem.hashtable += htBytes;
			
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

			for (size_t i = lo; i < hi; ++i) {
				u_kmer = kmers[i];
				kmer_t prefix = GET_PREFIX_SHIFTED(u_kmer);
				suffix_t suffix = GET_SUFFIX(u_kmer);
				
#ifdef USE_PREFETCH
				/*
				if (i + 2 * prefetch_dist < hi) {
					prefetch_kmer = kmers[i + 2 * prefetch_dist];
					kmer_t prefetch_prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
					
					const char* addr = static_cast<const char*>(&hashtables[prefetch_prefix]);
#ifdef WIN32
					_mm_prefetch(addr, _MM_HINT_T0);
#else
					__builtin_prefetch(data);
#endif
				}
				*/

				if (i + prefetch_dist < hi) {
					prefetch_kmer = kmers[i + prefetch_dist];
					kmer_t prefetch_prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
					suffix_t suffix = GET_SUFFIX(prefetch_kmer);
					hashtables[prefetch_prefix].prefetch(suffix);
				}
#endif
				// Check whether k-mer exists in a dictionary
				pattern_id_t* entry = hashtables[prefix].find(suffix);
				pattern_id_t p_id;

				if (entry == nullptr) {
					// do not add kmer to hashtable - just mark as to be added
					samplePatterns[to_add_id].first.kmer = u_kmer;
					--to_add_id;
				} else {
					p_id = *entry;

					samplePatterns[existing_id].first.pattern_id = p_id;
					samplePatterns[existing_id].second = entry;
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
void PrefixKmerDb::hashtableAdditionJob() {
	while (!this->queues.hashtableAddition.IsCompleted()) {
		HashtableAdditionTask task;

		if (this->queues.hashtableAddition.Pop(task)) {
			size_t lo = (*task.ranges)[task.block_id];
			size_t hi = (*task.ranges)[task.block_id + 1];

		
#ifdef USE_PREFETCH
			uint64_t prefetch_kmer;
#endif
			for (int j = lo; j < hi; ++j)
			{
#ifdef USE_PREFETCH
				if (j + prefetch_dist < hi) {
					prefetch_kmer = samplePatterns[kmers_to_add_to_HT[j + prefetch_dist]].first.kmer;

					kmer_t prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
					suffix_t suffix = GET_SUFFIX(prefetch_kmer);

					assert(prefix < hashtables.size());

					hashtables[prefix].prefetch(suffix);
				}
#endif
				int i = kmers_to_add_to_HT[j];
				kmer_t kmer = samplePatterns[i].first.kmer;

				kmer_t prefix = GET_PREFIX_SHIFTED(kmer);
				suffix_t suffix = GET_SUFFIX(kmer);

				auto i_kmer = hashtables[prefix].insert(suffix, 0); // First k-mer occurence - assign with temporary pattern 0
				samplePatterns[i].first.pattern_id = 0;
				samplePatterns[i].second = i_kmer;
			}

			semaphore.dec();

		}	
	}
}

// *****************************************************************************************
//
void PrefixKmerDb::hashtableJob() {

	while (!this->queues.hashtableAddition.IsCompleted()) {
		HashtableAdditionTask task;

		if (this->queues.hashtableAddition.Pop(task)) {
			// determine ranges
			const std::vector<kmer_t>& kmers = *(task.kmers);

			size_t lo = (*task.ranges)[task.block_id];
			size_t hi = (*task.ranges)[task.block_id + 1];

			kmer_t lo_prefix = GET_PREFIX_SHIFTED(kmers[lo]);
			kmer_t hi_prefix = GET_PREFIX_SHIFTED(kmers[hi - 1]);

			// generate histogram
			size_t begin = lo;
			kmer_t prevPrefix = lo_prefix;
			
			for (size_t i = lo + 1; i < hi; ++i) {
				kmer_t prefix = GET_PREFIX_SHIFTED(kmers[i]);
				if (prefix != prevPrefix) {
					prefixHistogram[prevPrefix] = i - begin;
					begin = i;
					prevPrefix = prefix;
				}
			}

			// add remanining
			prefixHistogram[hi_prefix] = hi - begin;
			
			// resize hastables
			size_t htBytes = 0;
			
			for (auto i = lo_prefix; i <= hi_prefix; ++i) {
				hashtables[i].reserve_for_additional(prefixHistogram[i]);
				htBytes += hashtables[i].get_bytes();
			}

		//	cout << "Block: " << task.block_id << ", lo_prefix: " << lo_prefix << ", hi_prefix: " << hi_prefix << endl << flush;


			// update memory statistics (atomic - no sync needed)
			mem.hashtable += htBytes;
			 
			size_t existing_id = lo;
			size_t to_add = 0;
			kmer_t u_kmer;
#ifdef USE_PREFETCH
			kmer_t prefetch_kmer;
			const size_t prefetch_dist = 48;
#endif

			for (size_t i = lo; i < hi; ++i) {
				u_kmer = kmers[i];
				kmer_t prefix = GET_PREFIX_SHIFTED(u_kmer);
				suffix_t suffix = GET_SUFFIX(u_kmer);

#ifdef USE_PREFETCH
				
				if (i + prefetch_dist < hi) {
					prefetch_kmer = kmers[i + prefetch_dist];
					kmer_t prefetch_prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
					suffix_t suffix = GET_SUFFIX(prefetch_kmer);
					hashtables[prefetch_prefix].prefetch(suffix);
				}
#endif
				// Check whether k-mer exists in a dictionary
				pattern_id_t* entry = hashtables[prefix].find(suffix);
				
				if (entry == nullptr) {
					entry = hashtables[prefix].insert(suffix, 0);
					++to_add;
				}
				
				samplePatterns[existing_id].first.pattern_id = *entry;
				samplePatterns[existing_id].second = entry;
				existing_id++;
			}

			kmersCount += to_add;

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
			int64_t deltaSize = 0; // get current pattern memory size (atomic)

			for (size_t i = lo; i < hi;) {
				size_t j;
				auto p_id = samplePatterns[i].first.pattern_id;

				// count k-mers from current sample with considered template 
				for (j = i + 1; j < hi; ++j) {
					if (p_id != samplePatterns[j].first.pattern_id) {
						break;
					}
				}
				size_t pid_count = j - i;

				if (patterns[p_id].get_num_kmers() == pid_count && !patterns[p_id].get_is_parrent()) {
					// Extend pattern - all k-mers with considered template exist in the analyzed sample
					deltaSize -= patterns[p_id].get_bytes();
					patterns[p_id].expand((sample_id_t)(task.sample_id));
					deltaSize += patterns[p_id].get_bytes();
				}
				else
				{
					// Generate new template 
					pattern_id_t local_pid = task.new_pid->fetch_add(1);

					threadPatterns[task.block_id].emplace_back(local_pid, pattern_t(patterns[p_id], p_id, task.sample_id, (uint32_t)pid_count));
					deltaSize += threadPatterns[task.block_id].back().second.get_bytes();

					if (p_id) {
						patterns[p_id].set_num_kmers(patterns[p_id].get_num_kmers() - pid_count);
					}

					for (size_t k = i; k < j; ++k) {
						*(samplePatterns[k].second) = local_pid;
					}

				}

				i = j;
			}

			// update memory statistics (atomic - no sync needed)
			mem.pattern += deltaSize;

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
	size_t n_kmers = kmers.size();
	
	//--------------------------------------------------------------------------
	// get prefix histogram (parallel)
	LOG_DEBUG << "Restructurizing hashtable (parallel)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	std::fill(prefixHistogram.begin(), prefixHistogram.end(), 0);
	samplePatterns.resize(n_kmers);

	mem.hashtable = 0;
	
	// prepare tasks
	size_t num_blocks = num_threads;
	size_t block = n_kmers / num_blocks;
	std::vector<size_t> ranges(num_blocks + 1, n_kmers);
	ranges[0] = 0;
	auto currentIndex = block;

	auto prefix_comparer = [this](kmer_t a, kmer_t b)->bool {
		return GET_PREFIX_SHIFTED(a) < GET_PREFIX_SHIFTED(b);
	};

	for (int tid = 0; tid < num_blocks - 1; ++tid) {
		auto it = std::upper_bound(
			kmers.begin() + currentIndex,
			kmers.end(),
			*(kmers.begin() + currentIndex - 1),
			prefix_comparer);

		size_t range = it - kmers.begin();

		ranges[tid + 1] = range;
		currentIndex = range + block;

		if (currentIndex >= kmers.size()) {
			break;
		}
	}

	semaphore.inc(num_blocks);

	for (int tid = 0; tid < num_blocks; ++tid) {
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.hashtableAddition.Push(HashtableAdditionTask{tid, &kmers, &ranges });
	}

	semaphore.waitForZero();
	times.hashtableResize += std::chrono::high_resolution_clock::now() - start;

#ifdef _DEBUG
	uint32_t histoSum = std::accumulate(prefixHistogram.begin(), prefixHistogram.end(), 0);
	if (histoSum != kmers.size()) {
		throw std::runtime_error("PrefixKmerDb::addKmers() - invalid histogram sum");
	}
#endif

	//--------------------------------------------------------------------------
	// sort in parallel
	LOG_DEBUG << "Sorting (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();

	ParallelSort(samplePatterns.data(), samplePatterns.size(), nullptr, 0, 0, num_threads);
	times.sort += std::chrono::high_resolution_clock::now() - start;
	
	//--------------------------------------------------------------------------
	// exdtend patterns
	LOG_DEBUG << "Extending patterns (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	std::atomic<size_t> new_pid(patterns.size());

	// calculate ranges
	std::fill(ranges.begin(), ranges.end(), n_kmers);
	ranges[0] = 0;
	block = n_kmers / num_threads;
	currentIndex = block;

	auto pid_comparer = [](const std::pair<kmer_or_pattern_t, pattern_id_t*>& a, const std::pair<kmer_or_pattern_t, pattern_id_t*>& b)->bool {
		return a.first.pattern_id < b.first.pattern_id;
	};

	for (int tid = 0; tid < num_threads - 1; ++tid) {
		auto it = std::upper_bound(
			samplePatterns.begin() + currentIndex,
			samplePatterns.end(),
			*(samplePatterns.begin() + currentIndex - 1),
			pid_comparer);

		size_t range = it - samplePatterns.begin();
		ranges[tid + 1] = range;
		currentIndex = range + block;

		if (currentIndex >= samplePatterns.size()) {
			break;
		}
	}

	// this should never happen
	if (ranges[num_threads] != n_kmers) {
		throw std::runtime_error("ERROR in FastKmerDb::addKmers(): Invalid ranges");
	}

	semaphore.inc(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.patternExtension.Push(PatternExtensionTask{ tid, sampleId, &ranges, &new_pid, nullptr });
	}

	// wait for the task to complete
	semaphore.waitForZero();
	times.extension += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// patterns insertion
	LOG_DEBUG << "Inserting patterns (serial)..." << endl;
	
	// extend by 1.5 on reallocation
	if (patterns.capacity() < new_pid) {
		patterns.reserve(new_pid * 3 / 2);
	}

	patterns.resize(new_pid);

	for (int tid = 0; tid < num_threads; ++tid) {
		for (auto& tp : threadPatterns[tid]) {
			patterns[tp.first] = std::move(tp.second);
		}
	}

	return sampleId;
}

/*
sample_id_t PrefixKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers, kmerLength, fraction);
	size_t n_kmers = kmers.size();

	//--------------------------------------------------------------------------
	// get prefix histogram (parallel)
	LOG_DEBUG << "Restructurizing hashtable (parallel)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	std::fill(prefixHistogram.begin(), prefixHistogram.end(), 0);
	mem.hashtable = 0;
	// prepare tasks
	semaphore.inc(num_threads);
	internalSempahores[0].inc(num_threads);
	internalSempahores[1].inc(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
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

	//--------------------------------------------------------------------------
	// find for kmers in parallel
	std::vector<size_t> num_existing_kmers(num_threads);
	samplePatterns.resize(n_kmers);

	LOG_DEBUG << "Finding kmers (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	// prepare tasks
	semaphore.inc(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.hashtableSearch.Push(DictionarySearchTask{ tid, num_threads, &kmers, &num_existing_kmers[tid] });
	}
	// wait for the task to complete
	semaphore.waitForZero();
	times.hashtableFind += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// add kmers to hashtable sequentially
	LOG_DEBUG << "Adding kmers (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();

	kmers_to_add_to_HT.clear();
	for (int tid = 0; tid < num_threads; ++tid) {
		size_t n_kmers = kmers.size();
		size_t block = n_kmers / num_threads;
		size_t lo = tid * block;
		size_t hi = (tid == num_threads - 1) ? n_kmers : lo + block;

		for (size_t i = num_existing_kmers[tid]; i < hi; ++i)
			kmers_to_add_to_HT.push_back(i);
	}

	size_t n_kmers_to_add = kmers_to_add_to_HT.size();
	kmersCount += n_kmers_to_add;

	// calculate ranges
	std::vector<size_t> ranges(num_threads + 1, n_kmers_to_add);
	ranges[0] = 0;
	size_t block = n_kmers_to_add / num_threads;
	auto currentIndex = block;

	auto prefix_comparer = [this](int a, int b)->bool {
		return GET_PREFIX_SHIFTED(this->samplePatterns[a].first.kmer) < GET_PREFIX_SHIFTED(this->samplePatterns[b].first.kmer);
	};

	for (int tid = 0; tid < num_threads - 1; ++tid) {
		auto it = std::upper_bound(
			kmers_to_add_to_HT.begin() + currentIndex,
			kmers_to_add_to_HT.end(),
			*(kmers_to_add_to_HT.begin() + currentIndex - 1),
			prefix_comparer);

		size_t range = it - kmers_to_add_to_HT.begin();

		ranges[tid + 1] = range;
		currentIndex = range + block;

		if (currentIndex >= kmers_to_add_to_HT.size()) {
			break;
		}
	}


	semaphore.inc(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.hashtableAddition.Push(HashtableAdditionTask{ tid, &ranges });
	}

	semaphore.waitForZero();

	/*
	#ifdef USE_PREFETCH
	uint64_t prefetch_kmer;
	#endif
	for (int j = 0; j < n_kmers_to_add; ++j)
	{
	#ifdef USE_PREFETCH
	if (j + prefetch_dist < n_kmers_to_add) {
	prefetch_kmer = samplePatterns[kmers_to_add_to_HT[j + prefetch_dist]].first.kmer;

	kmer_t prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
	suffix_t suffix = GET_SUFFIX(prefetch_kmer);

	assert(prefix < hashtables.size());

	hashtables[prefix].prefetch(suffix);
	}
	#endif
	int i = kmers_to_add_to_HT[j];
	kmer_t kmer = samplePatterns[i].first.kmer;

	kmer_t prefix = GET_PREFIX_SHIFTED(kmer);
	suffix_t suffix = GET_SUFFIX(kmer);

	auto i_kmer = hashtables[prefix].insert(suffix, 0); // First k-mer occurence - assign with temporary pattern 0
	samplePatterns[i].first.pattern_id = 0;
	samplePatterns[i].second = i_kmer;
	}

	
	times.hashtableAdd += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// sort in parallel
	LOG_DEBUG << "Sorting (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();

	ParallelSort(samplePatterns.data(), samplePatterns.size(), nullptr, 0, 0, num_threads);
	times.sort += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// exdtend patterns
	LOG_DEBUG << "Extending patterns (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	std::atomic<size_t> new_pid(patterns.size());

	// calculate ranges
	std::fill(ranges.begin(), ranges.end(), n_kmers);
	ranges[0] = 0;
	block = n_kmers / num_threads;
	currentIndex = block;

	auto pid_comparer = [](const std::pair<kmer_or_pattern_t, pattern_id_t*>& a, const std::pair<kmer_or_pattern_t, pattern_id_t*>& b)->bool {
		return a.first.pattern_id < b.first.pattern_id;
	};

	for (int tid = 0; tid < num_threads - 1; ++tid) {
		auto it = std::upper_bound(
			samplePatterns.begin() + currentIndex,
			samplePatterns.end(),
			*(samplePatterns.begin() + currentIndex - 1),
			pid_comparer);

		size_t range = it - samplePatterns.begin();
		ranges[tid + 1] = range;
		currentIndex = range + block;

		if (currentIndex >= samplePatterns.size()) {
			break;
		}
	}

	// this should never happen
	if (ranges[num_threads] != n_kmers) {
		throw std::runtime_error("ERROR in FastKmerDb::addKmers(): Invalid ranges");
	}

	semaphore.inc(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		queues.patternExtension.Push(PatternExtensionTask{ tid, sampleId, &ranges, &new_pid, nullptr });
	}

	// wait for the task to complete
	semaphore.waitForZero();
	times.extension += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// patterns insertion
	LOG_DEBUG << "Inserting patterns (serial)..." << endl;

	// extend by 1.5 on reallocation
	if (patterns.capacity() < new_pid) {
		patterns.reserve(new_pid * 3 / 2);
	}

	patterns.resize(new_pid);

	for (int tid = 0; tid < num_threads; ++tid) {
		for (auto& tp : threadPatterns[tid]) {
			patterns[tp.first] = std::move(tp.second);
		}
	}

	return sampleId;
}

*/

// *****************************************************************************************
//
void PrefixKmerDb::serialize(std::ofstream& file) const {

	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<suffix_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
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
				file.write(buffer, bufpos * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));
				bufpos = 0;
			}
		}
		// write remaining ht elements
		file.write(reinterpret_cast<const char*>(&bufpos), sizeof(size_t));
		file.write(buffer, bufpos * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));
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

	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<suffix_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());

	LOG_VERBOSE << "Loading general info..." << endl;

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

	LOG_VERBOSE << "Loading kmer hashtables..." << endl;

	// load number of hashmaps
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	hashtables.resize(temp);

	// load all hashtables
	for (int i = 0; i < hashtables.size(); ++i) {
		cout << i << ",";
		auto& ht = hashtables[i];
		// loal ht size
		file.read(reinterpret_cast<char*>(&temp), sizeof(temp));

		// load ht elements
		ht.clear();
		ht.reserve_for_additional(temp);

		size_t readCount = 0;
		while (readCount < temp) {
			size_t portion = 0;
			file.read(reinterpret_cast<char*>(&portion), sizeof(size_t));
			file.read(buffer, portion * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));

			for (size_t j = 0; j < portion; ++j) {
				ht.insert(hashtableBuffer[j].key, hashtableBuffer[j].val);
			}
			readCount += portion;
		}
	}

	cout << endl;

	if (!file) {
		return false;
	}

	cout << "Loading patterns..." << endl;

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

