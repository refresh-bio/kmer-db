/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <algorithm>
#include <cstdio>
#include <string>
#include <map>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <chrono>
#include <stack>
#include <thread>
#include <mutex>
#include <atomic>
#include <iterator>
#include <numeric>
#include <emmintrin.h>
#include <immintrin.h>
#include "instrset.h"

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "log.h"

#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif
#include <omp.h>


#define USE_PREFETCH
#define ALL_STATS

using namespace std;

const size_t FastKmerDb::ioBufferBytes = (2 << 29); //512MB buffer 
//const size_t FastKmerDb::ioBufferBytes = 16000000; //16MB buffer 
//const size_t FastKmerDb::ioBufferBytes = 100000; //100KB buffer 


// *****************************************************************************************
//
FastKmerDb::FastKmerDb(int _num_threads, size_t cacheBufferMb) :
	kmers2patternIds(), 
//	repeatedKmers(),
	dictionarySearchQueue(1), 
	patternExtensionQueue(1),
	num_threads(_num_threads > 0 ? _num_threads : std::thread::hardware_concurrency()),
	cacheBufferMb(cacheBufferMb) {

	avx2_present = instrset_detect() >= 8;

	patternBytes = 0;
	patterns.reserve(2 << 20);
	patterns.push_back(pattern_t());
	
//	threadRepeatedKmers.resize(num_threads);
	threadPatterns.resize(num_threads);
	
	for (int tid = 0; tid < num_threads; ++tid) {
//		threadRepeatedKmers[tid].reserve(2 << 20);
		threadPatterns[tid].reserve(2 << 20);
	}

	dictionarySearchWorkers.resize(num_threads);
	patternExtensionWorkers.resize(num_threads);
	
	for (auto& t : dictionarySearchWorkers) {
		t = std::thread([this]() {
			while (!this->dictionarySearchQueue.IsCompleted()) {
				DictionarySearchTask task;

				if (this->dictionarySearchQueue.Pop(task)) {
					LOG_DEBUG << "Block " << task.block_id << " started" << endl;
					const std::vector<kmer_t>& kmers = *(task.kmers);
					
					size_t n_kmers = kmers.size();
					size_t block_size = n_kmers / task.num_blocks;
					size_t lo = task.block_id * block_size;
					size_t hi = (task.block_id == task.num_blocks - 1) ? n_kmers : lo + block_size;

//					threadRepeatedKmers[task.block_id].clear();
//					threadRepeatedKmers[task.block_id].reserve(hi - lo);

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
						if (u_kmer & KMER_MSB) {
							u_kmer &= ~KMER_MSB; // disable non-uniqueness flag
//							threadRepeatedKmers[task.block_id].push_back(u_kmer); // mark as sample-repeated
						} 

#ifdef USE_PREFETCH
						if (i + prefetch_dist < hi)
						{
							prefetch_kmer = kmers[i + prefetch_dist] & (~KMER_MSB);
							kmers2patternIds.prefetch(prefetch_kmer);
						}
#endif
						// Check whether k-mer exists in a dictionary
						auto i_kmer = kmers2patternIds.find(u_kmer);
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
		});
	}
	

	for (auto& t : patternExtensionWorkers) {
		t = std::thread([this]() {
			while (!this->patternExtensionQueue.IsCompleted()) {
				PatternExtensionTask task;

				if (this->patternExtensionQueue.Pop(task)) {
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
							patterns[p_id].expand((sample_id_t) (task.sample_id));
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
		});
	}

}

// *****************************************************************************************
//
FastKmerDb::~FastKmerDb() {
	dictionarySearchQueue.MarkCompleted();
	for (auto& t : dictionarySearchWorkers) {
		t.join();
	}

	patternExtensionQueue.MarkCompleted();
	for (auto& t : patternExtensionWorkers) {
		t.join();
	}
}

// *****************************************************************************************
//
sample_id_t FastKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction)
{
#ifdef USE_PREFETCH
	const size_t prefetch_dist = 48;
#endif
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers, kmerLength, fraction);

	size_t n_kmers = kmers.size();
	
	LOG_DEBUG << "Restructurizing hashtable (serial)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();

	// to prevent hashtable restructuring during k-mer addition
	kmers2patternIds.reserve_for_additional(n_kmers);
	// repeatedKmers.reserve_for_additional(n_kmers);
	
	samplePatterns.resize(n_kmers);

	std::vector<size_t> num_existing_kmers(num_threads);
	hashtableResizeTime += std::chrono::high_resolution_clock::now() - start;
	
	// find for kmers in parallel
	LOG_DEBUG << "Finding kmers (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> threads(num_threads);
	// prepare tasks
	for (int tid = 0; tid < num_threads; ++tid) {
		semaphore.inc();
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		dictionarySearchQueue.Push(DictionarySearchTask{ tid, num_threads, &kmers, &num_existing_kmers[tid] });
	}
	// wait for the task to complete
	semaphore.waitForZero();
	hashtableFindTime += std::chrono::high_resolution_clock::now() - start;
	
	// add kmers to hashtable sequentially
	LOG_DEBUG << "Adding kmers (serial)..." << endl;
	start = std::chrono::high_resolution_clock::now();

	kmers_to_add_to_HT.clear();
	for (int tid = 0; tid < num_threads; ++tid) {
		size_t n_kmers = kmers.size();
		size_t block = n_kmers / num_threads;
		size_t lo = tid * block;
		size_t hi = (tid == num_threads - 1) ? n_kmers : lo + block;

		// copy repeated k-mers from thread to global collection
/*		for (auto it = threadRepeatedKmers[tid].cbegin(); it != threadRepeatedKmers[tid].cend(); ++it) {
			if (!repeatedKmers.is_free(*it)) {
				repeatedKmers.insert(*it);
			}
		}
*/		
		for (size_t i = num_existing_kmers[tid]; i < hi; ++i)
			kmers_to_add_to_HT.push_back(i);
	}

	size_t n_kmers_to_add = kmers_to_add_to_HT.size();
#ifdef USE_PREFETCH
	uint64_t prefetch_kmer;
#endif
	for (int j = 0; j < n_kmers_to_add; ++j)
	{
#ifdef USE_PREFETCH
		if (j + prefetch_dist < n_kmers_to_add) {
			prefetch_kmer = samplePatterns[kmers_to_add_to_HT[j + prefetch_dist]].first;
			kmers2patternIds.prefetch(prefetch_kmer);
		}
#endif
		int i = kmers_to_add_to_HT[j];
		auto i_kmer = kmers2patternIds.insert(samplePatterns[i].first, 0); // First k-mer occurence - assign with temporary pattern 0
		samplePatterns[i].first = 0;
		samplePatterns[i].second = i_kmer;
	}

	hashtableAddTime += std::chrono::high_resolution_clock::now() - start;
	
	LOG_DEBUG << "Sorting (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();

	ParallelSort(samplePatterns.data(), samplePatterns.size(), nullptr, 0, 0, num_threads);

	
	sortTime += std::chrono::high_resolution_clock::now() - start;
	LOG_VERBOSE << "sort time: " << sortTime.count() << "  " << samplePatterns.size() << endl;


	LOG_DEBUG << "Extending patterns (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	std::atomic<size_t> new_pid(patterns.size());

	// calculate ranges
	std::vector<size_t> ranges(num_threads + 1, n_kmers);
	ranges[0] = 0;
	size_t block = n_kmers / num_threads;

	auto currentIndex = block;

	std::vector<size_t> threadBytes(num_threads, 0);
	auto pid_comparer = [](const std::pair<pattern_id_t, pattern_id_t*>& a, const std::pair<pattern_id_t, pattern_id_t*>& b)->bool {
		return a.first < b.first;
	};

	for (int tid = 0; tid < num_threads-1; ++tid) {
		auto it = std::upper_bound(
			samplePatterns.begin() + currentIndex,
			samplePatterns.end(), 
			*(samplePatterns.begin() + currentIndex - 1), 
			pid_comparer);

		size_t range = it - samplePatterns.begin();
		
	
		ranges[tid + 1] = range;
		currentIndex =  range + block;

		if (currentIndex >= samplePatterns.size()) {
			break;
		}
	}
	
	// this should never happen
	if (ranges[num_threads] != n_kmers) {
		throw std::runtime_error("ERROR in FastKmerDb::addKmers(): Invalid ranges");
	}


	for (int tid = 0; tid < num_threads; ++tid) {
		semaphore.inc();
		LOG_DEBUG << "Block " << tid << " scheduled" << endl;
		patternExtensionQueue.Push(PatternExtensionTask{ tid, sampleId, &ranges, &new_pid, &threadBytes });
	}
		
	// wait for the task to complete
	semaphore.waitForZero();
	extensionTime += std::chrono::high_resolution_clock::now() - start;

	LOG_DEBUG << "Inserting kmers (serial)..." << endl;
	for (int tid = 0; tid < threads.size(); ++tid) {
		patternBytes += threadBytes[tid];
	}

	// extend by 1.5 on reallocation
	if (patterns.capacity() < new_pid) {
		patterns.reserve(new_pid * 3 / 2);
	}

	patterns.resize(new_pid);

	for (int tid = 0; tid < threads.size(); ++tid) {
		for (auto& tp : threadPatterns[tid]) {
			patterns[tp.first] = std::move(tp.second);
		}
	}

	return sampleId;
}

// *****************************************************************************************
//
void FastKmerDb::mapKmers2Samples(uint64_t kmer, std::vector<sample_id_t>& samples) const {
	
	// find corresponding pattern id
	auto p_id = kmers2patternIds.find(kmer);
	const auto& p = patterns[*p_id];
	samples.resize(p.get_num_samples());
	
	p.decodeSamples(samples.data());

}

// *****************************************************************************************
//
void FastKmerDb::serialize(std::ofstream& file) const {
	
	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<kmer_t, pattern_id_t>::item_t> hastableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hastableBuffer.data());

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

	// store hashmap
	temp = kmers2patternIds.get_size();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));
	
	// write ht elements in portions
	size_t bufpos = 0;
	for (auto it = kmers2patternIds.cbegin(); it < kmers2patternIds.cend(); ++it) {
		if (kmers2patternIds.is_free(*it)) {
			continue;
		}

		hastableBuffer[bufpos++] = *it;
		if (bufpos == numHastableElements) {
			file.write(reinterpret_cast<const char*>(&bufpos), sizeof(size_t));
			file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
			bufpos = 0;
		}
	}
	// write remaining ht elements
	file.write(reinterpret_cast<const char*>(&bufpos), sizeof(size_t));
	file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
	
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
bool FastKmerDb::deserialize(std::ifstream& file) {
	
	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<kmer_t, pattern_id_t>::item_t> hastableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hastableBuffer.data());

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

	// load hashtable
	temp = kmers2patternIds.get_size();
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));

	kmers2patternIds.clear();
	kmers2patternIds.reserve_for_additional(temp);

	size_t readCount = 0;
	while (readCount < temp) {
		size_t portion = 0;
		file.read(reinterpret_cast<char*>(&portion), sizeof(size_t));
		file.read(buffer, portion * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));

		for (size_t j = 0; j < portion; ++j) {
			kmers2patternIds.insert(hastableBuffer[j].key, hastableBuffer[j].val);
		}
		readCount += portion;
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

// *****************************************************************************************
//
void FastKmerDb::savePatterns(std::ofstream& file) const {
	
	std::vector<uint32_t> aux(getSamplesCount());
	
	for (int i = 0; i < patterns.size(); ++i) {
		const auto& p = patterns[i];
		file << i << ": " << p.get_parent_id() << " | ";
		p.decodeSamples(aux.data());
		std::copy(aux.begin(), aux.begin() + p.get_num_local_samples(), std::ostream_iterator<uint32_t>(file, " "));
		file << endl;
	}

}

// *****************************************************************************************
//
void FastKmerDb::calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix) //const 
{
	int samples_count = getSamplesCount();
	matrix.resize(samples_count);
	matrix.clear();
	
	size_t bufsize = cacheBufferMb * 1000000 / sizeof(uint32_t);
	std::vector<uint32_t> patternsBuffer(bufsize);
	uint32_t* currentPtr;

	std::vector<uint32_t*> rawPatterns(bufsize);
	std::vector<std::tuple<sample_id_t, uint32_t, uint32_t>> sample2pattern(bufsize);

	std::vector<std::thread> workers_matrix(num_threads);
	std::vector<std::thread> workers_decomp(num_threads);
	std::vector<std::thread> workers_sample2patterns(num_threads);
	std::vector<std::thread> workers_histogram(num_threads);

	int first_pid;

	// Set number of k-mers to internal nodes
	int num_patterns = patterns.size();
	for (int i = num_patterns-1; i > 0; --i)
	{
		int parent_id = patterns[i].get_parent_id();

		if (parent_id >= 0)
			patterns[parent_id].add_num_kmers(patterns[i].get_num_kmers());
	}

	// determine ranges of blocks processed by threads 
	int no_ranges = num_threads * 16;
	std::vector<size_t> workerRanges(no_ranges + 1, sample2pattern.size());
	std::vector<int> v_range_ids;

	for (int i = no_ranges - 1; i >= 0; --i)
		v_range_ids.push_back(i);

	CRegisteringQueue<int> tasks_matrix_queue(1);
	Semaphore semaphore_matrix;
	std::atomic<uint64_t> numAdditions(0);

	// ranges of patterns to decompress
	std::vector<int> v_range_patterns(no_ranges + 1);
	CRegisteringQueue<pair<int, uint32_t>> tasks_decomp_queue(1);
	Semaphore semaphore_decomp;

	CRegisteringQueue<int> tasks_sample2patterns_queue(1);
	Semaphore semaphore_sample2patterns;

	Semaphore semaphore_hist_first;		// semaphore for first stage in histogram calculation
	Semaphore semaphore_hist_second;		// semaphore for second stage in histogram calculation

	int no_hist_parts = num_threads * 1;
	std::vector<std::vector<int>> hist_sample_ids(no_hist_parts, std::vector<int>(samples_count, 0));
	CRegisteringQueue<int> tasks_histogram_queue(1);
	std::vector<int> hist_boundary_values(num_threads, 0);
	std::vector<pair<int, uint32_t>> v_tmp;
	std::vector<int> v_tmp_int;
	std::vector<int> v_hist_ids(num_threads);
	for (int i = 0; i < num_threads; ++i)
		v_hist_ids[i] = i;

	// Decompress patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_decomp[tid] = std::thread([&sample2pattern, &v_range_patterns, &tasks_decomp_queue, &semaphore_decomp, &rawPatterns, &first_pid, this, &patternsBuffer, &hist_sample_ids]
		{
			while (!tasks_decomp_queue.IsCompleted())
			{
				pair<int, uint32_t> decomp_task;

				if (tasks_decomp_queue.Pop(decomp_task))
				{
					int part_id = decomp_task.first;
					auto &my_hist_sample_ids = hist_sample_ids[part_id];
					my_hist_sample_ids.clear();
					my_hist_sample_ids.resize(this->getSamplesCount(), 0);

					int f_pid = v_range_patterns[part_id];
					int l_pid = v_range_patterns[part_id + 1];
					auto currentPtr = patternsBuffer.data() + decomp_task.second;

					for (int pid = f_pid; pid < l_pid; ++pid)
					{
						const auto& pattern = this->patterns[pid];

						currentPtr += pattern.get_num_samples();
						uint32_t* out = currentPtr;		// start from the end

						// decode all samples from pattern and its parents
						int64_t current_id = pid;
						while (current_id >= 0) {
							const auto& cur = patterns[current_id];
							auto parent_id = cur.get_parent_id();

							if (parent_id >= 0)
								_mm_prefetch((const char*)(patterns.data() + parent_id), _MM_HINT_T0);

							out -= cur.get_num_local_samples();
							cur.decodeSamples(out);

							current_id = parent_id;
						}
						rawPatterns[pid - first_pid] = out; // begin of unpacked pattern

						int num_samples = pattern.get_num_samples();
						int num_local_samples = pattern.get_num_local_samples();
						for (int i = num_samples - num_local_samples; i < num_samples; ++i)
							++my_hist_sample_ids[out[i]];
					}

					semaphore_decomp.dec();
				}
			}
		});
	}
	
	// Sample2patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_sample2patterns[tid] = std::thread([&sample2pattern, &v_range_patterns, &tasks_sample2patterns_queue, &semaphore_sample2patterns, &rawPatterns, &first_pid, this, &patternsBuffer, &hist_sample_ids]
		{
			while (!tasks_sample2patterns_queue.IsCompleted())
			{
				int part_id;

				if (tasks_sample2patterns_queue.Pop(part_id))
				{
					auto &my_hist_sample_ids = hist_sample_ids[part_id];

					int f_pid = v_range_patterns[part_id];
					int l_pid = v_range_patterns[part_id + 1];

					for (int pid = f_pid; pid < l_pid; ++pid)
					{
						const auto& pattern = this->patterns[pid];

						auto out = rawPatterns[pid - first_pid];
						int num_samples = pattern.get_num_samples();
						int num_local_samples = pattern.get_num_local_samples();

						for (int i = num_samples - num_local_samples; i < num_samples; ++i)
						{
							int pos = my_hist_sample_ids[out[i]]++;

							sample2pattern[pos] = make_tuple(out[i], pid - first_pid, i);
						}
					}

					semaphore_sample2patterns.dec();
				}
			}
		});
	}

	// increment array elements in threads
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_matrix[tid] = std::thread([&sample2pattern, &rawPatterns, &workerRanges, &matrix, this, &tasks_matrix_queue, &first_pid, &semaphore_matrix, &numAdditions] {
			uint64_t localAdditions = 0;

			int range_id;

			while (!tasks_matrix_queue.IsCompleted())
			{
				if (tasks_matrix_queue.Pop(range_id))
				{
					for (int id = workerRanges[range_id]; id < workerRanges[range_id + 1]; ) {
						int Si = get<0>(sample2pattern[id]);
						uint32_t *row = matrix[Si];
						while (id < workerRanges[range_id + 1] && get<0>(sample2pattern[id]) == Si) {
							int local_pid = get<1>(sample2pattern[id]);
							const auto& pattern = patterns[local_pid + first_pid];

							uint32_t* rawData = rawPatterns[local_pid];
							int num_samples = get<2>(sample2pattern[id]);
							uint32_t to_add = pattern.get_num_kmers();

#ifdef ALL_STATS
							localAdditions += num_samples;
#endif

							auto *p = rawData;

							row_add(row, p, num_samples, to_add, avx2_present);	
							++id;
						}
					}
					semaphore_matrix.dec();
				}
			}
			numAdditions.fetch_add(localAdditions);
		});
	}

	// calculate histogram
	for (int tid = 0; tid < num_threads; ++tid)
		workers_histogram[tid] = std::thread([&hist_boundary_values, &tasks_histogram_queue, samples_count, this, &v_tmp, 
			&hist_sample_ids, &semaphore_hist_first, &semaphore_hist_second] {
			int task_id;

			while (!tasks_histogram_queue.IsCompleted())
			{
				if (tasks_histogram_queue.Pop(task_id))
				{
					int first_sample = samples_count * task_id / num_threads;
					int last_sample = samples_count * (task_id + 1) / num_threads;

					int sum = 0;
					int max_j = v_tmp.size();
					for (int i = first_sample; i < last_sample; ++i)
						for (int j = 0; j < max_j; ++j)
						{
							int x = hist_sample_ids[j][i];
							hist_sample_ids[j][i] = sum;
							sum += x;
						}

					hist_boundary_values[task_id] = sum;
					semaphore_hist_first.dec_notify_all();
					semaphore_hist_first.waitForZero();

					int to_add = 0;
					for (int i = 0; i < task_id; ++i)
						to_add += hist_boundary_values[i];

					for (int j = 0; j < max_j; ++j)
						for (int i = first_sample; i < last_sample; ++i)
							hist_sample_ids[j][i] += to_add;
					semaphore_hist_second.dec();
				}
			}
		});

	// process all patterns in blocks determined by buffer size
	std::cout << std::endl;

	double decomp_time = 0;
	double hist_time = 0;
	double patterns_time = 0;

	int pid_to_cout = 0;
	for (int pid = 0; pid < patterns.size(); ) {
		if (pid >= pid_to_cout)
		{
			std::cout << pid << " of " << patterns.size() << "\r";
			fflush(stdout);
			pid_to_cout += 50000;
		}

		first_pid = pid;
		currentPtr = patternsBuffer.data();
		size_t samplesCount = 0;

		auto t1 = std::chrono::high_resolution_clock::now();
		
		int part_size = bufsize / no_hist_parts + 1;
		int next_boundary = 0;
		int part_id = 0;

		v_tmp.clear();
		v_tmp_int.clear();

		while (pid < patterns.size() && samplesCount + patterns[pid].get_num_samples() < bufsize && part_id < no_hist_parts) {
			const auto& pattern = patterns[pid];
			
			if (samplesCount >= next_boundary)
			{
				next_boundary += part_size;
				v_range_patterns[part_id] = pid;
				v_tmp.push_back(make_pair(part_id, (uint32_t) samplesCount));
				v_tmp_int.push_back(part_id);

				++part_id;
			}

			samplesCount += pattern.get_num_samples();
			++pid;
		}

		v_range_patterns[part_id] = pid;
		semaphore_decomp.inc(v_tmp.size());
		tasks_decomp_queue.PushRange(v_tmp);
		semaphore_decomp.waitForZero();

		int last_pid = pid;
		int num_samples = getSamplesCount();
		
		auto t2 = std::chrono::high_resolution_clock::now();

		semaphore_hist_first.inc(num_threads);
		semaphore_hist_second.inc(num_threads);
		tasks_histogram_queue.PushRange(v_hist_ids);
		semaphore_hist_second.waitForZero();

		auto t3 = std::chrono::high_resolution_clock::now();

		// generate sample to pattern mapping
		int num_items = std::accumulate(hist_boundary_values.begin(), hist_boundary_values.end(), 0);
		sample2pattern.resize(num_items);		// resize to number of items in the current fragment

		semaphore_sample2patterns.inc(v_tmp_int.size());
		tasks_sample2patterns_queue.PushRange(v_tmp_int);
		semaphore_sample2patterns.waitForZero();

		auto t4 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> dt1 = t2 - t1;
		std::chrono::duration<double> dt2 = t3 - t2;
		std::chrono::duration<double> dt3 = t4 - t3;
		
		decomp_time += dt1.count();
		hist_time += dt2.count();
		patterns_time += dt3.count();

		// determine ranges of blocks processed by threads 
		workerRanges.assign(no_ranges + 1, sample2pattern.size());
		workerRanges[0] = 0;
		size_t workerBlock = sample2pattern.size() / no_ranges;

		auto currentIndex = workerBlock;

		for (int rid = 0; rid < no_ranges - 1; ++rid) {
			// make sure no single row of matrix is updated by multiple workers
			auto it = std::upper_bound(
				sample2pattern.begin() + currentIndex,
				sample2pattern.end(),
				*(sample2pattern.begin() + currentIndex - 1),
				[](auto x, auto y) {return get<0>(x) < get<0>(y); });	// Necessary, because we are looking for end of sample data

			size_t range = it - sample2pattern.begin();
			workerRanges[rid + 1] = range;
			currentIndex = range + workerBlock;

			if (currentIndex >= sample2pattern.size())
				break;
		}

		// this should never happen
		if (workerRanges[no_ranges] != sample2pattern.size()) {
			throw std::runtime_error("ERROR in FastKmerDb::calculateSimilarity(): Invalid ranges");
		}

		semaphore_matrix.inc(no_ranges);
		tasks_matrix_queue.PushRange(v_range_ids);
		semaphore_matrix.waitForZero();
	}

	tasks_matrix_queue.MarkCompleted();
	tasks_decomp_queue.MarkCompleted();
	tasks_sample2patterns_queue.MarkCompleted();
	tasks_histogram_queue.MarkCompleted();

	for (auto & w : workers_matrix)
		w.join();

	for (auto & w : workers_decomp)
		w.join();

	for (auto & w : workers_sample2patterns)
		w.join();

	for (auto &w : workers_histogram)
		w.join();

	cout << "Decomp   time: " << decomp_time << endl;
	cout << "Hist     time: " << hist_time << endl;
	cout << "Sample2p time: " << patterns_time << endl;

#ifdef ALL_STATS
	cout << "Number of additions:" << numAdditions << endl;
#endif
}

// *****************************************************************************************
//
void  FastKmerDb::calculateSimilarity(const std::vector<kmer_t>& kmers, std::vector<uint32_t>& similarities) const {
	similarities.resize(this->getSamplesCount(), 0);

	std::vector<std::vector<uint32_t>> localSimilarities(num_threads, std::vector<uint32_t>(this->getSamplesCount()));

	std::unordered_map<pattern_id_t, int32_t> patterns2count;

	std::chrono::duration<double> dt;
	auto start = std::chrono::high_resolution_clock::now();

	// iterate over kmers in analyzed sample
	for (const auto& kmer: kmers) {
		// check if kmer exists in a database
		auto entry = kmers2patternIds.find(kmer);

		if (entry != nullptr) {
			auto pid = *entry;
			const auto& pattern = patterns[pid];

			if (pattern.get_num_kmers() == 0)
				continue;

			++patterns2count[pid];
		}
	}

	std::vector<std::pair<pattern_id_t, int32_t>> patterns2countVector(patterns2count.size());
	int i = 0;
	uint64_t tot_pattern_lengths = 0;
	for (const auto& entry : patterns2count) {
		tot_pattern_lengths += patterns[entry.first].get_num_samples();
	}

	vector<uint32_t> range_boundaries(num_threads + 1, patterns2count.size());
	uint64_t sum_pattern_lengths = 0;
	uint64_t part_length = (tot_pattern_lengths + num_threads - 1) / num_threads;
	uint32_t range_id = 1;
	for (auto &entry : patterns2count)
	{
		sum_pattern_lengths += patterns[entry.first].get_num_samples();
		patterns2countVector[i++] = entry;

		if (sum_pattern_lengths > range_id * part_length)
			range_boundaries[range_id++] = i;
	}
	range_boundaries[0] = 0;

	patterns2count.clear();


	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_VERBOSE << "Pattern listing time: " << dt.count() << endl;

	start = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> workers(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, tid, &patterns2countVector, &localSimilarities, &range_boundaries]() {
			std::vector<uint32_t> samples(this->getSamplesCount());
			
			size_t n_patterns = patterns2countVector.size();
			size_t block_size = n_patterns / num_threads;
			size_t lo = range_boundaries[tid];
			size_t hi = range_boundaries[tid+1];
			auto &my_localSimilarities = localSimilarities[tid];

			for (int id = lo; id < hi; ++id) {
				if (id + 1 < hi)
					_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);

				auto pid = patterns2countVector[id].first;
				const auto& pattern = patterns[pid];
				int num_samples = pattern.get_num_samples();
				int to_add = patterns2countVector[id].second;

				uint32_t* out = samples.data() + pattern.get_num_samples(); // start from the end

				int64_t current_id = pid;
				while (current_id >= 0) {
					const auto& cur = patterns[current_id];

					out -= cur.get_num_local_samples();
					cur.decodeSamples(out);

					current_id = cur.get_parent_id();
				}

				auto *p = samples.data();

				int i;
				for (i = 0; i + 4 <= num_samples; i += 4)
				{
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
				}
				num_samples -= i;

				switch (num_samples)
				{
				case 3: my_localSimilarities[*p++] += to_add;
				case 2: my_localSimilarities[*p++] += to_add;
				case 1: my_localSimilarities[*p++] += to_add;
				}

			}
		});
	}

	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid].join();
		std::transform(localSimilarities[tid].begin(), localSimilarities[tid].end(), similarities.begin(), similarities.begin(), [](uint32_t a, uint32_t b)->uint32_t{ 
			return a + b;
		});
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_VERBOSE << "Pattern unpacking time: " << dt.count() << endl;
}
