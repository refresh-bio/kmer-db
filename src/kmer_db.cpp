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
#include <cassert>

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
FastKmerDb::FastKmerDb(int _num_threads) :
//	repeatedKmers(),
	dictionarySearchQueue(1), 
	patternExtensionQueue(1),
	num_threads(_num_threads) {

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
	for (int tid = 0; tid < num_threads; ++tid) {
		patternBytes += threadBytes[tid];
	}

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

// *****************************************************************************************
//
void FastKmerDb::mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const {
	
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

