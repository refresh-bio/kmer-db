#include "prefix_kmer_db.h"
#include "prefix_kmer_db.h"

#include "../libs/refresh/sort/lib/pdqsort_par.h"

#include "log.h"

#include <algorithm>

#include <numeric>
#include <cassert>

#define USE_PREFETCH
#define ALL_STATS

using namespace std;

// *****************************************************************************************
//
PrefixKmerDb::PrefixKmerDb(int _num_threads) 
	: AbstractKmerDb(), num_threads(_num_threads)
{
	patterns.reserve(2 << 20);
	patterns.push_back(pattern_t());

	threadPatterns.resize(num_threads);
	for (auto& tp : threadPatterns) { 
		tp.reserve(2 << 10); 
	}
}

// *****************************************************************************************
//
PrefixKmerDb::~PrefixKmerDb() {
	queues.hashtableAddition.AddTasks(std::vector<HashtableTask>());
	for (auto& t : workers.hashtableAddition) {
		t.join();
	}

	queues.patternExtension.AddTasks(std::vector<PatternTask>());
	for (auto& t : workers.patternExtension) {
		t.join();
	}
}

// *****************************************************************************************
//
void PrefixKmerDb::initialize(uint32_t kmerLength, double fraction, AlphabetType alphabetType) {
	AbstractKmerDb::initialize(kmerLength, fraction, alphabetType);

	// temporary instance
	std::shared_ptr<Alphabet> alphabet(AlphabetFactory::instance().create(alphabetType));
	
	int prefixBits = (int)kmerLength * alphabet->bitsPerSymbol - SUFFIX_BITS;
	
	if (prefixBits < 8) {
		prefixBits = 8;
	}

	size_t binsCount = 1ULL << prefixBits;

	hashtables.resize(binsCount);
}

// *****************************************************************************************
//
void PrefixKmerDb::hashtableJobATP() {

	HashtableTask task;

	size_t local_hashtableBytes = 0;

	while (this->queues.hashtableAdditionATP.Pop(task)) 
	{
		// determine ranges
		const kmer_t* kmers = task.kmers;

		uint32_t lo = task.lo;
		uint32_t hi = task.hi;

		if (lo != hi) {

			kmer_t lo_prefix = GET_PREFIX_SHIFTED(kmers[lo]);
			kmer_t hi_prefix = GET_PREFIX_SHIFTED(kmers[hi - 1]);

#ifdef COLLECT_DETAILED_TIMES
			auto start = std::chrono::high_resolution_clock::now();
#endif
			// generate histogram
			uint32_t begin = lo;
			kmer_t prevPrefix = lo_prefix;

			// resize hastables
			size_t htBytes = 0;
			int reallocsCount = 0;

			uint32_t dif = 0;

			for (uint32_t i = lo + 1; i < hi; ++i) {
				kmer_t prefix = GET_PREFIX_SHIFTED(kmers[i]);
				if (prefix != prevPrefix) {
//					prefixHistogram[prevPrefix] = i - begin;
					dif = i - begin;

					// ---
//					reallocsCount += hashtables[prevPrefix].reserve_for_additional(prefixHistogram[prevPrefix]);
					reallocsCount += hashtables[prevPrefix].reserve_for_additional(dif);
					htBytes += hashtables[prevPrefix].get_bytes();
					// ---

					begin = i;
					prevPrefix = prefix;
				}
			}

			// add remanining
			dif  = hi - begin;

			reallocsCount += hashtables[hi_prefix].reserve_for_additional(dif);
			htBytes += hashtables[hi_prefix].get_bytes();

			local_hashtableBytes += htBytes;

#ifdef COLLECT_DETAILED_TIMES
			times.hashtableResize_worker += std::chrono::high_resolution_clock::now() - start;
#endif

			//	LOG_DEBUG << "Block: " << task.block_id << ", lo_prefix: " << lo_prefix << ", hi_prefix: " << hi_prefix << endl ;

#ifdef COLLECT_DETAILED_TIMES
			start = std::chrono::high_resolution_clock::now();
#endif

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

				if (entry->val == ht.empty_value) {
					ht.insert(suffix, 0, entry);
					++to_add;
				}

				samplePatterns[existing_id].first.pattern_id = entry->val;
				samplePatterns[existing_id].second = &(entry->val);
				existing_id++;
			}

			kmersCount += to_add;
#ifdef COLLECT_DETAILED_TIMES
			times.hashtableFind_worker += std::chrono::high_resolution_clock::now() - start;
#endif
		}
	}

	stats.hashtableBytes += local_hashtableBytes;
}

// *****************************************************************************************
//
void PrefixKmerDb::patternJobATP() {
	PatternTask task;

	size_t local_patternBytes = 0;

	while (this->queues.patternExtensionATP.Pop(task)) 
	{
		uint32_t lo = task.lo;
		uint32_t hi = task.hi;
			
	//	LOG_DEBUG << "Pattern job " << task.block_id << " started (" << lo << "-" << hi << ")" << endl ;
			
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

		local_patternBytes += deltaSize;
	}

	stats.patternBytes += local_patternBytes;
}

// *****************************************************************************************
//
sample_id_t PrefixKmerDb::addKmers(
	const std::string& sampleName,
	const kmer_t* kmers,
	uint32_t kmersCount,
	uint32_t kmerLength,
	double fraction,
	AlphabetType alphabetType,
	refresh::active_thread_pool& atp)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers, kmersCount, kmerLength, fraction, alphabetType, atp);
	uint32_t n_kmers = static_cast<uint32_t>(kmersCount);
	
	if (n_kmers == 0)
	{
		LOG_NORMAL("Empty sample: " << sampleName << endl);
		return sampleId;
	}

	//--------------------------------------------------------------------------
	// get prefix histogram (parallel)
	
	LOG_DEBUG("Hashtable resizing, searching, and adding (parallel)..." << endl);
#ifdef COLLECT_DETAILED_TIMES
	auto start = std::chrono::high_resolution_clock::now();
#endif

	samplePatterns.resize(n_kmers);

	stats.hashtableBytes = 0;
	
	// prepare tasks
	uint32_t num_blocks = (uint32_t)num_threads * 8;
//	uint32_t block = std::max(n_kmers / num_blocks, 1u);
	uint32_t block = std::max(n_kmers / num_blocks, 64u);

//	bool use_busy_wait = num_threads >= 8;

	// prepare tasks
	std::vector<HashtableTask> hashtableTasks(num_blocks, HashtableTask{0, 0, n_kmers, kmers, n_kmers });
	uint32_t currentHi = 0;

	auto prefix_comparer = [this](kmer_t a, kmer_t b)->bool {
		return GET_PREFIX_SHIFTED(a) < GET_PREFIX_SHIFTED(b);
	};

	uint32_t tid = 0;
	for (tid = 0; tid < num_blocks && currentHi < n_kmers; ++tid) {
		// set block id and low bound
		hashtableTasks[tid].block_id = tid;
		hashtableTasks[tid].lo = (tid == 0) ? 0 : hashtableTasks[tid - 1].hi;
	
		currentHi = hashtableTasks[tid].lo + block;

		// check if it makes sense to search for upper bound
		if (currentHi < n_kmers && tid < num_blocks - 1) {
			kmer_t ref = *(kmers + currentHi - 1);
			auto it = std::upper_bound(kmers + currentHi, kmers + n_kmers, ref, prefix_comparer);

			hashtableTasks[tid].hi = (uint32_t)(it - kmers); // this is always positive
		}	
	}

	hashtableTasks.resize(tid);

	// Queue works like a stack, so large items should be at the end
/*	std::sort(hashtableTasks.begin(), hashtableTasks.end(), [](const HashtableTask& x, const HashtableTask& y)->bool {
		return (x.hi - x.lo) < (y.hi - y.lo);
	});*/

	refresh::sort::pdqsort(hashtableTasks.begin(), hashtableTasks.end(), [](const HashtableTask& x, const HashtableTask& y)->bool {
//	stable_sort(hashtableTasks.begin(), hashtableTasks.end(), [](const HashtableTask& x, const HashtableTask& y)->bool {
		return (x.hi - x.lo) < (y.hi - y.lo);
		});

	stats.hashtableJobsImbalance += (double)(hashtableTasks.front().hi - hashtableTasks.front().lo) * num_blocks / n_kmers;

	refresh::active_thread_pool_state hashtableAddition_state;

	queues.hashtableAdditionATP.Push(hashtableTasks);

	int n_ht_jobs = std::min<int>(num_threads - 1, (int) hashtableTasks.size() - 1);

	for (int i = 0; i < n_ht_jobs; ++i)
		atp.launch([&, this] {hashtableJobATP(); }, &hashtableAddition_state);
	hashtableJobATP();
	if(n_ht_jobs)
		hashtableAddition_state.busy_wait();

//	semaphore.waitForZero();
#ifdef COLLECT_DETAILED_TIMES
	times.hashtableProcess += std::chrono::high_resolution_clock::now() - start;
#endif
	
	//--------------------------------------------------------------------------
	// sort in parallel
	LOG_DEBUG("Sorting (parallel)..." << endl);
#ifdef COLLECT_DETAILED_TIMES
	start = std::chrono::high_resolution_clock::now();
#endif

	refresh::sort::pdqsort_branchless_tp(refresh::sort::pdqsort_adjust_threads(samplePatterns.size(), num_threads), samplePatterns.begin(), samplePatterns.end(), 
		[](const auto& a, const auto& b) {return a.first.pattern_id < b.first.pattern_id; }, atp);

#ifdef COLLECT_DETAILED_TIMES
	times.sort += std::chrono::high_resolution_clock::now() - start;
#endif
	
	//--------------------------------------------------------------------------
	// exdtend patterns
	LOG_DEBUG("Extending patterns (parallel)..." << endl);
#ifdef COLLECT_DETAILED_TIMES
	start = std::chrono::high_resolution_clock::now();
#endif
	std::atomic<pattern_id_t> new_pid((pattern_id_t)patterns.size());

	// prepare tasks
	num_blocks = num_threads;
	block = std::max(n_kmers / num_blocks, 1u);
	std::vector<PatternTask> patternTasks(num_blocks, PatternTask{ 0, n_kmers, n_kmers, sampleId, &new_pid});
	patternTasks[0].lo = 0;

	auto pid_comparer = [](const std::pair<kmer_or_pattern_t, pattern_id_t*>& a, const std::pair<kmer_or_pattern_t, pattern_id_t*>& b)->bool {
		return a.first.pattern_id < b.first.pattern_id;
	};

	auto currentIndex = block;
	for (uint32_t tid = 0; tid < num_blocks - 1; ++tid) {
		auto it = std::upper_bound(
			samplePatterns.begin() + currentIndex,
			samplePatterns.end(),
			*(samplePatterns.begin() + currentIndex - 1),
			pid_comparer);

		uint32_t range = (uint32_t)(it - samplePatterns.begin());

		patternTasks[tid + 1].lo = patternTasks[tid].hi = range;
		patternTasks[tid + 1].block_id = tid + 1;
		
		block = std::max((n_kmers - range) / (num_blocks - tid - 1), 1u);
		
		currentIndex = range + block;

		if (currentIndex >= samplePatterns.size()) {
			break;
		}
	}

	patternTasks.erase(
		std::find_if(patternTasks.begin(), patternTasks.end(), [](const PatternTask& t)->bool { return t.hi == t.lo;  }),
		patternTasks.end());

//	std::sort(patternTasks.begin(), patternTasks.end(), [](const auto& x, const auto& y) {return x.hi - x.lo < y.hi - y.lo; });
	refresh::sort::pdqsort(patternTasks.begin(), patternTasks.end(), [](const auto& x, const auto& y) {return x.hi - x.lo < y.hi - y.lo; });
//	stable_sort(patternTasks.begin(), patternTasks.end(), [](const auto& x, const auto& y) {return x.hi - x.lo < y.hi - y.lo; });

	queues.patternExtensionATP.Push(patternTasks);

	refresh::active_thread_pool_state patternExtension_state;
	
	int n_pat_jobs = std::min<int>(num_threads - 1, (int) patternTasks.size() - 1);

	for (int i = 0; i < n_pat_jobs; ++i)
		atp.launch([&, this] {patternJobATP(); }, &patternExtension_state);
	patternJobATP(); 
	
	if (n_pat_jobs)
		patternExtension_state.busy_wait();

#ifdef COLLECT_DETAILED_TIMES
	times.extension += std::chrono::high_resolution_clock::now() - start;
#endif

	//--------------------------------------------------------------------------
	// patterns insertion
	LOG_DEBUG("Moving patterns to global collection (serial)..." << endl);
	
	// extend by 1.5 on reallocation
	if (patterns.capacity() < (size_t)new_pid) {
		patterns.reserve(new_pid * 3 / 2);
	}

	patterns.resize(new_pid);

	for (size_t tid = 0; tid < patternTasks.size(); ++tid) {
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
	save(file, formatWord);
	
	// save kmer length and fraction
	save(file, kmerLength);
	save(file, fraction);
	save(file, startFraction);
	save(file, alphabetType);
	save(file, isInitialized);
	save(file, kmersCount);

	// store number of samples
	size_t temp = getSamplesCount();
	save(file, temp);

	// store sample info 
	for (size_t i = 0; i < sampleNames.size(); ++i) {
		temp = sampleKmersCount[i]; // store kmer count
		save(file, temp);

		const string& s = sampleNames[i]; // store name
		temp = s.size();
		save(file, temp);
		file.write(s.data(), temp);
	}

	// store number of hashmaps
	LOG_NORMAL("Storing k-mer hashtables (" << (rawHashtables ? "raw" : "compressed"s) << ")..." << endl);
	auto time = std::chrono::high_resolution_clock::now();

	temp = hashtables.size();
	save(file, temp);

	size_t percent = 0;

	// store all hashmaps
	for (size_t i = 0; i < hashtables.size(); ++i) {
		
		if (i * 100 / hashtables.size() > percent) {
			LOG_NORMAL("\r" << percent << "%");
			++percent;
		}

		auto& ht = hashtables[i];

		if (rawHashtables) {
			ht.serialize(file, hashtableBuffer.data(), hashtableBuffer.size());
		}
		else {
			// store ht size
			temp = ht.get_size();
			save(file, temp);

			if (temp > 0) {
				// write ht elements in portions
				size_t accum = 0;
				size_t bufpos = 0;
				for (auto it = ht.cbegin(); it < ht.cend(); ++it) {
					if (ht.is_free(*it)) {
						continue;
					}

					hashtableBuffer[bufpos++] = *it;
					if (bufpos == numHastableElements) {
						save(file, bufpos);
						file.write(buffer, bufpos * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));
						accum += bufpos;
						bufpos = 0;
					}
				}

				// write remaining ht elements
				save(file, bufpos);
				file.write(buffer, bufpos * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));
				accum += bufpos;
				
				if (accum != ht.get_size()) {
					throw std::runtime_error("Assertion error : HT size does not match the number of elements");
				}
			}
		}
	}
	std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - time;
	LOG_NORMAL("\r" << hashtables.size() << "/" << hashtables.size() << " hashtables stored in " << dt.count() << " s");
	LOG_NORMAL(endl);

	// write patterns in portions
	LOG_NORMAL("Storing patterns..." << endl);

	time = std::chrono::high_resolution_clock::now();
	temp = patterns.size();
	save(file, temp);

	char * currentPtr = buffer;

	percent = 0;

	for (size_t pid = 0; pid < patterns.size(); ++pid) {
		
		if (pid * 100 / patterns.size() > percent) {
			LOG_NORMAL("\r" << percent << "%");
			++percent;
		}
		
		if (currentPtr + patterns[pid].get_bytes() > buffer + IO_BUFFER_BYTES) {
			size_t blockSize = currentPtr - buffer;
			save(file, blockSize); // write size of block to facilitate deserialization
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
	save(file, blockSize); // write size of block to facilitate deserialization
	file.write(buffer, blockSize);

	dt = std::chrono::high_resolution_clock::now() - time;
	LOG_NORMAL("\r" << patterns.size() << "/" << patterns.size() << " patterns stored in " << dt.count() << " s");
	LOG_NORMAL(endl);
}

// *****************************************************************************************
//
bool PrefixKmerDb::deserialize(std::ifstream& file, DeserializationMode mode) {

	size_t numHastableElements = IO_BUFFER_BYTES / sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<suffix_t, pattern_id_t>::item_t> hashtableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hashtableBuffer.data());

	// load format word
	bool rawHashtables = false;

	uint64_t formatWord = 0;
	load(file, formatWord);
	if (formatWord & SERIALIZATION_RAW_HASHTABLES) {
		rawHashtables = true;
	}

	// load kmer length and fraction
	load(file, kmerLength);
	load(file, fraction);
	load(file, startFraction);
	load(file, alphabetType);
	load(file, isInitialized);
	load(file, kmersCount);

	// load sample info
	size_t temp;
	load(file, temp);
	sampleNames.resize(temp);
	sampleKmersCount.resize(temp);

	for (size_t i = 0; i < sampleNames.size(); ++i) {
		load(file, temp);  // load kmer count
		sampleKmersCount[i] = temp;

		string& s = sampleNames[i];
		load(file, temp);  // load sample name
		file.read(buffer, temp);
		s.assign(buffer, temp); 
	}

	if (!file) {
		return false;
	}

	if (mode == DeserializationMode::SamplesOnly) {
		return true;
	}

	if (mode != DeserializationMode::SkipHashtables) {
		LOG_NORMAL("Loading k-mer hashtables (" << (rawHashtables ? "raw" : "compressed"s) << ")..." << endl);
	}
	
	auto time = std::chrono::high_resolution_clock::now();

	// load number of hashmaps
	load(file, temp);
	suffix_kmers.resize(temp);
	hashtables.resize(temp);			

	size_t percent = 0;

	// load all hashtables
	for (size_t i = 0; i < hashtables.size(); ++i) {
		
		if (mode != DeserializationMode::SkipHashtables && (i * 100 / hashtables.size() > percent)) {
			LOG_NORMAL("\r" << percent << "%");
			++percent;
		}

		auto& ht = hashtables[i];

		if (rawHashtables) {

			if (mode == DeserializationMode::CompactedHashtables) {
				ht.deserialize_into_vector(file, hashtableBuffer.data(), hashtableBuffer.size(), suffix_kmers[i], false);
			}
			else {
				ht.deserialize(file, hashtableBuffer.data(), hashtableBuffer.size(), (mode == DeserializationMode::SkipHashtables));
			}
		}
		else {
			// load ht size
			size_t temp;
			load(file, temp);

			if (temp > 0) {

				// load ht elements
				if (mode == DeserializationMode::CompactedHashtables) {
					suffix_kmers[i].reserve(temp);
				}
				else {
					ht.clear();
					ht.reserve_for_additional(temp);
				}

				size_t readCount = 0;
				while (readCount < temp) {
					size_t portion = 0;
					load(file, portion);

					if (mode == DeserializationMode::SkipHashtables) {
						file.seekg(portion, ios_base::cur);
					}
					else {
						file.read(reinterpret_cast<char*>(buffer), portion * sizeof(hash_map_lp<suffix_t, pattern_id_t>::item_t));

						if (mode == DeserializationMode::CompactedHashtables) {
							for (size_t j = 0; j < portion; ++j) {
								suffix_kmers[i].emplace_back(hashtableBuffer[j].key, hashtableBuffer[j].val);
							}
						}
						else {
							for (size_t j = 0; j < portion; ++j) {
								ht.insert(hashtableBuffer[j].key, hashtableBuffer[j].val);
							}
						}
					}
					readCount += portion;
				}
			}
		}
	}
	
	std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - time;
	if (mode != DeserializationMode::SkipHashtables) {
		LOG_NORMAL("\r" << hashtables.size() << "/" << hashtables.size() << " hashtables loaded in " << dt.count() << " s");
		LOG_NORMAL(endl);
	}

	if (!file) {
		return false;
	}

	LOG_NORMAL("Loading patterns..." << endl);
	time = std::chrono::high_resolution_clock::now();

	// load patterns
	load(file, temp);
	patterns.clear();
	patterns.resize(temp);

	percent = 0;

	size_t pid = 0;
	while (pid < patterns.size()) {
		size_t blockSize;
		load(file, blockSize);
		file.read(buffer, blockSize);

		char * currentPtr = buffer;
		while (currentPtr < buffer + blockSize) {

			if (pid * 100 / patterns.size() > percent) {
				LOG_NORMAL("\r" << percent << "%");
				++percent;
			}

			currentPtr = patterns[pid].unpack(currentPtr);
			++pid;
		}
	}
	dt = std::chrono::high_resolution_clock::now() - time;
	LOG_NORMAL("\r" << patterns.size() << "/" << patterns.size() << " patterns loaded in " << dt.count() << " s");
	LOG_NORMAL(endl);

	if (!file) {
		return false;
	}

	return true;
}

// *****************************************************************************************
//
void PrefixKmerDb::savePatterns(std::ofstream& file) const {

	std::vector<uint32_t> aux(getSamplesCount());

	for (size_t i = 0; i < patterns.size(); ++i) {
		const auto& p = patterns[i];
		file << i << ": " << p.get_parent_id() << " | ";
		p.decodeSamples(aux.data());
		std::copy(aux.begin(), aux.begin() + p.get_num_local_samples(), std::ostream_iterator<uint32_t>(file, " "));
		file << endl;
	}
}
