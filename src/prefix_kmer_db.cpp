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
	num_threads(_num_threads)
{
	patterns.reserve(2 << 20);
	patterns.push_back(pattern_t());

	threadPatterns.resize(num_threads);
	for (auto& tp : threadPatterns) { 
		tp.reserve(2 << 10); 
	}
	
	workers.hashtableAddition.resize(num_threads);
	workers.patternExtension.resize(num_threads);

	for (auto& t : workers.hashtableAddition) {
		t = std::thread(&PrefixKmerDb::hashtableJob, this);
	}

	for (auto& t : workers.patternExtension) {
		t = std::thread(&PrefixKmerDb::patternJob, this);
	}

	//stats.workerAdditions.resize(num_threads);
	//stats.workerReallocs.resize(num_threads);
}

// *****************************************************************************************
//
PrefixKmerDb::~PrefixKmerDb() {
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
void PrefixKmerDb::hashtableJob() {

	while (!this->queues.hashtableAddition.IsCompleted()) {
		HashtableTask task;

		if (this->queues.hashtableAddition.Pop(task)) {
			// determine ranges
			const kmer_t* kmers = task.kmers;
		
			uint32_t lo = task.lo;
			uint32_t hi = task.hi;

			if (lo != hi) {

				kmer_t lo_prefix = GET_PREFIX_SHIFTED(kmers[lo]);
				kmer_t hi_prefix = GET_PREFIX_SHIFTED(kmers[hi - 1]);

				auto start = std::chrono::high_resolution_clock::now();
				// generate histogram
				size_t begin = lo;
				kmer_t prevPrefix = lo_prefix;

				for (uint32_t i = lo + 1; i < hi; ++i) {
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
				int reallocsCount = 0;
				for (auto i = lo_prefix; i <= hi_prefix; ++i) {
					reallocsCount += hashtables[i].reserve_for_additional(prefixHistogram[i]);
					htBytes += hashtables[i].get_bytes();
				}

				// update memory statistics (atomic - no sync needed)
				stats.hashtableBytes += htBytes;
			//	stats.workerReallocs[task.block_id] += reallocsCount;

				times.hashtableResize_worker += std::chrono::high_resolution_clock::now() - start;

				LOG_DEBUG << "Block: " << task.block_id << ", lo_prefix: " << lo_prefix << ", hi_prefix: " << hi_prefix << endl << flush;

				start = std::chrono::high_resolution_clock::now();

				uint32_t existing_id = lo;
				uint32_t to_add = 0;
				kmer_t u_kmer;
#ifdef USE_PREFETCH
				kmer_t prefetch_kmer;
#endif

				for (uint32_t i = lo; i < hi; ++i) {
					u_kmer = kmers[i];
					kmer_t prefix = GET_PREFIX_SHIFTED(u_kmer);
					suffix_t suffix = GET_SUFFIX(u_kmer);

#ifdef USE_PREFETCH

					if (i + PREFETCH_DIST < hi) {
						prefetch_kmer = kmers[i + PREFETCH_DIST];
						kmer_t prefetch_prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
						suffix_t suffix = GET_SUFFIX(prefetch_kmer);
						hashtables[prefetch_prefix].prefetch(suffix);
					}
#endif
					// Check whether k-mer exists in a dictionary
					auto& ht = hashtables[prefix];
					auto* entry = ht.find_item(suffix);

					if (entry->key == ht.empty_key) {
						ht.insert(suffix, 0, entry);
						++to_add;
					}

					samplePatterns[existing_id].first.pattern_id = entry->val;
					samplePatterns[existing_id].second = &(entry->val);
					existing_id++;
				}

				kmersCount += to_add;
				times.hashtableFind_worker += std::chrono::high_resolution_clock::now() - start;
			}

			this->semaphore.dec();
		}
	}

}

// *****************************************************************************************
//
void PrefixKmerDb::patternJob() {
	while (!this->queues.patternExtension.IsCompleted()) {
		PatternTask task;

		if (this->queues.patternExtension.Pop(task)) {
			LOG_DEBUG << "Pattern job " << task.block_id << " started" << endl;
			
			uint32_t lo = task.lo;
			uint32_t hi = task.hi;
			sample_id_t sampleId = (sample_id_t)(task.sample_id);
			
			threadPatterns[task.block_id].clear();
			int64_t deltaSize = 0; // get current pattern memory size (atomic)

			for (uint32_t i = lo; i < hi;) {
				uint32_t j;
				auto p_id = samplePatterns[i].first.pattern_id;

				// count k-mers from current sample with considered template 
				for (j = i + 1; j < hi; ++j) {
					if (p_id != samplePatterns[j].first.pattern_id) {
						break;
					}
				}
				uint32_t pid_count = j - i;

				if (patterns[p_id].get_num_kmers() == pid_count && !patterns[p_id].get_is_parrent()) {
					// Extend pattern - all k-mers with considered template exist in the analyzed sample
					deltaSize -= patterns[p_id].get_bytes();
					patterns[p_id].expand(sampleId);
					deltaSize += patterns[p_id].get_bytes();
				}
				else
				{
					// Generate new template 
					pattern_id_t local_pid = task.new_pid->fetch_add(1);

					threadPatterns[task.block_id].emplace_back( local_pid, pattern_t(patterns[p_id], p_id, sampleId, pid_count) );
					deltaSize += threadPatterns[task.block_id].back().second.get_bytes();

					if (p_id) {
						patterns[p_id].set_num_kmers(patterns[p_id].get_num_kmers() - pid_count);
					}

					for (uint32_t k = i; k < j; ++k) {
						*(samplePatterns[k].second) = local_pid;
					}
				}

				i = j;
			}

			// update memory statistics (atomic - no sync needed)
			stats.patternBytes += deltaSize;

			LOG_DEBUG << "Pattern job " << task.block_id << " finished" << endl;
			this->semaphore.dec();
		}
	}
}


// *****************************************************************************************
//
sample_id_t PrefixKmerDb::addKmers(
	const std::string& sampleName,
	const kmer_t* kmers,
	size_t kmersCount,
	uint32_t kmerLength,
	double fraction)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers, kmersCount, kmerLength, fraction);
	uint32_t n_kmers = static_cast<uint32_t>(kmersCount);
	
//	std::fill(stats.workerReallocs.begin(), stats.workerReallocs.end(), 0);
//	std::fill(stats.workerAdditions.begin(), stats.workerAdditions.end(), 0);

	//--------------------------------------------------------------------------
	// get prefix histogram (parallel)
	LOG_DEBUG << "Hashtable resizing, searching, and adding (parallel)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	std::fill(prefixHistogram.begin(), prefixHistogram.end(), 0);
	samplePatterns.resize(n_kmers);

	stats.hashtableBytes = 0;
	
	// prepare tasks
	size_t num_blocks = num_threads * 8;
	size_t block = n_kmers / num_blocks;

	// prepare tasks
	std::vector<HashtableTask> hashtableTasks(num_blocks, HashtableTask{0, 0, n_kmers, kmers, n_kmers });
	auto currentHi = 0;

	auto prefix_comparer = [this](kmer_t a, kmer_t b)->bool {
		return GET_PREFIX_SHIFTED(a) < GET_PREFIX_SHIFTED(b);
	};

	int tid = 0;
	for (tid = 0; tid < num_blocks && currentHi < n_kmers; ++tid) {
		// set block id and low bound
		hashtableTasks[tid].block_id = tid;
		hashtableTasks[tid].lo = (tid == 0) ? 0 : hashtableTasks[tid - 1].hi;

		currentHi = hashtableTasks[tid].lo + block;

		// check if it makes sense to search for upper bound
		if (currentHi < n_kmers) {

			auto it = std::upper_bound(kmers + currentHi, kmers + n_kmers,
				*(kmers + currentHi - 1), prefix_comparer);

			hashtableTasks[tid].hi = it - kmers;
		}	
	}

	hashtableTasks.resize(tid);

	std::sort(hashtableTasks.begin(), hashtableTasks.end(), [](const HashtableTask& x, const HashtableTask& y)->bool {
		return (x.hi - x.lo) > (y.hi - y.lo);
	});

	stats.hashtableJobsImbalance += (double)(hashtableTasks.front().hi - hashtableTasks.front().lo) * num_blocks / n_kmers;

	semaphore.inc(hashtableTasks.size());
	for (int tid = 0; tid < hashtableTasks.size(); ++tid) {
		LOG_DEBUG << "Hashtable job " << tid << " scheduled" << endl;
		queues.hashtableAddition.Push(hashtableTasks[tid]);
	}

	semaphore.waitForZero();
	times.hashtableProcess += std::chrono::high_resolution_clock::now() - start;
	
#ifdef _DEBUG
	uint32_t histoSum = std::accumulate(prefixHistogram.begin(), prefixHistogram.end(), 0);
	if (histoSum != n_kmers) {
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
	std::atomic<pattern_id_t> new_pid(patterns.size());

	// prepare tasks
	num_blocks = num_threads;
	block = n_kmers / num_blocks;
	std::vector<PatternTask> patternTasks(num_blocks, PatternTask{ 0, n_kmers, n_kmers, sampleId, &new_pid});
	patternTasks[0].lo = 0;

	auto pid_comparer = [](const std::pair<kmer_or_pattern_t, pattern_id_t*>& a, const std::pair<kmer_or_pattern_t, pattern_id_t*>& b)->bool {
		return a.first.pattern_id < b.first.pattern_id;
	};

	/*
	currentHi = 0;

	tid = 0;
	for (tid = 0; tid < num_blocks && currentHi < n_kmers; ++tid) {
		// set block id and low bound
		patternTasks[tid].block_id = tid;
		patternTasks[tid].lo = (tid == 0) ? 0 : patternTasks[tid - 1].hi;

		currentHi = patternTasks[tid].lo + block;

		// check if it makes sense to search for upper bound
		if (currentHi < n_kmers) {

			auto it = std::upper_bound(samplePatterns.begin() + currentHi, samplePatterns.end(),
				*(samplePatterns.begin() + currentHi - 1), pid_comparer);

			patternTasks[tid].hi = it - samplePatterns.begin();
		}
	}

	patternTasks.resize(tid);

	*/

	auto currentIndex = block;
	for (int tid = 0; tid < num_blocks - 1; ++tid) {
		auto it = std::upper_bound(
			samplePatterns.begin() + currentIndex,
			samplePatterns.end(),
			*(samplePatterns.begin() + currentIndex - 1),
			pid_comparer);

		size_t range = it - samplePatterns.begin();

		patternTasks[tid + 1].lo = patternTasks[tid].hi = range;
		patternTasks[tid + 1].block_id = tid + 1;
		currentIndex = range + block;

		if (currentIndex >= samplePatterns.size()) {
			break;
		}
	}

	patternTasks.erase(
		std::find_if(patternTasks.begin(), patternTasks.end(), [](const PatternTask& t)->bool { return t.hi == t.lo;  }),
		patternTasks.end());

	semaphore.inc(patternTasks.size());
	for (int tid = 0; tid < patternTasks.size(); ++tid) {
		LOG_DEBUG << "Pattern job " << tid << " scheduled" << endl;
		queues.patternExtension.Push(patternTasks[tid]);
	}

	// wait for the task to complete
	semaphore.waitForZero();
	times.extension += std::chrono::high_resolution_clock::now() - start;

	//--------------------------------------------------------------------------
	// patterns insertion
	LOG_DEBUG << "Moving patterns to global collection (serial)..." << endl;
	
	// extend by 1.5 on reallocation
	if (patterns.capacity() < new_pid) {
		patterns.reserve(new_pid * 3 / 2);
	}

	patterns.resize(new_pid);

	for (int tid = 0; tid < patternTasks.size(); ++tid) {
		for (auto& tp : threadPatterns[tid]) {
			patterns[tp.first] = std::move(tp.second);
		}
	}

	return sampleId;
}

// *****************************************************************************************
//
void PrefixKmerDb::serialize(std::ofstream& file, bool rawHashtables) const {

	size_t numHastableElements = IO_BUFFER_BYTES / sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<suffix_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());
	
	// generate format word
	uint64_t formatWord = 0;
	if (rawHashtables) {
		formatWord |= SERIALIZATION_RAW_HASHTABLES;
	}
	file.write(reinterpret_cast<const char*>(&formatWord), sizeof(formatWord));
	
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

		if (rawHashtables) {
			// store ht capacity
			temp = ht.get_capacity();
			file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

			// write ht elements in portions
			size_t bufpos = 0;
			for (;;) {
				
	
			}
		}
		else {
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
	}

	// write patterns in portions
	temp = patterns.size();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

	char * currentPtr = buffer;
	for (int pid = 0; pid < patterns.size(); ++pid) {
		if (currentPtr + patterns[pid].get_bytes() > buffer + IO_BUFFER_BYTES) {
			size_t blockSize = currentPtr - buffer;
			file.write(reinterpret_cast<const char*>(&blockSize), sizeof(size_t)); // write size of block to facilitate deserialization
			file.write(buffer, blockSize);
			currentPtr = buffer;
		}

		currentPtr = patterns[pid].pack(currentPtr);
		// this should never happen
		if (currentPtr > buffer + IO_BUFFER_BYTES) {
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

	size_t numHastableElements = IO_BUFFER_BYTES / sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<suffix_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());

	cout << "Loading general info..." << endl;

	// load format word
	bool rawHashatable = false;

	uint64_t formatWord = 0;
	file.read(reinterpret_cast<char*>(&formatWord), sizeof(formatWord));
	if (formatWord & SERIALIZATION_RAW_HASHTABLES) {
		rawHashatable = true;
	}


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

	cout << "Loading kmer hashtables..." << endl;

	// load number of hashmaps
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	hashtables.resize(temp);

	// load all hashtables
	for (int i = 0; i < hashtables.size(); ++i) {
		if ((i + 1) % 10 == 0) {
			cout << "\r" << i + 1 << "/" << hashtables.size() << "...                      " << std::flush;
		}
		auto& ht = hashtables[i];
		
		if (rawHashatable) {


		}
		else {
			// load ht size
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

			if ((pid + 1) % 1000 == 0) {
				cout << "\r" << pid + 1 << "/" << patterns.size() << "...                      " << std::flush;
			}

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


// *****************************************************************************************
//
void PrefixKmerDb::savePatterns(std::ofstream& file) const {

	std::vector<uint32_t> aux(getSamplesCount());

	for (int i = 0; i < patterns.size(); ++i) {
		const auto& p = patterns[i];
		file << i << ": " << p.get_parent_id() << " | ";
		p.decodeSamples(aux.data());
		std::copy(aux.begin(), aux.begin() + p.get_num_local_samples(), std::ostream_iterator<uint32_t>(file, " "));
		file << endl;
	}

}