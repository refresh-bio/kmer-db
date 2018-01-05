// kmer_db.cpp : Defines the entry point for the console application.
//

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

/****************************************************************************************************************************************************/



/****************************************************************************************************************************************************/

sample_id_t NaiveKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers);
	
	uint64_t u_kmer;
#ifdef USE_PREFETCH
	uint64_t prefetch_kmer;
	const size_t prefetch_dist = 48;
#endif
	auto n_kmers = kmers.size();
	
	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	kmers2patternIds.reserve_for_additional(n_kmers);
	patterns.reserve(patterns.size() + n_kmers);

	for (size_t i = 0; i < n_kmers; ++i) {
		u_kmer = kmers[i];
#ifdef USE_PREFETCH
		if (i + prefetch_dist < n_kmers) {
			prefetch_kmer = kmers[i + prefetch_dist];
			kmers2patternIds.prefetch(prefetch_kmer);
		}
#endif

		// ***** Sprawdzanie w slowniku czy taki k-mer juz istnieje
		auto patternId = kmers2patternIds.find(u_kmer);

		if (patternId == nullptr) {
			patterns.push_back(std::vector<sample_id_t>(1, sampleId));
			kmers2patternIds.insert(u_kmer, patterns.size() - 1);
		}
		else {
			patterns[*patternId].push_back(sampleId);
		}

		patternBytes += sizeof(sample_id_t);
	}

	return sampleId;
}


void NaiveKmerDb::mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) const {
	samples.clear();
	auto v = kmers2patternIds.find(kmer);

	if (v != nullptr) {
		samples = patterns[*v];
	}
}


void NaiveKmerDb::calculateSimilarity(Array<uint32_t>& matrix) const {
	matrix.resize(getSamplesCount(), getSamplesCount());
	matrix.clear();

	for (auto it = kmers2patternIds.cbegin(); it < kmers2patternIds.cend(); ++it) {
		
		if (kmers2patternIds.is_free(*it)) {
			continue;
		}

		auto & samples = patterns[it->val];
		
		for (int i = 0; i < samples.size() - 1; ++i) {
			auto& Si = samples[i];
			for (int j = i + 1; j < samples.size(); ++j) {
				auto& Sj = samples[j];
				++matrix[Si][Sj];
				++matrix[Sj][Si];
			}
		}
	}
}

std::map<std::vector<sample_id_t>, size_t> NaiveKmerDb::getPatternsStatistics() const {
	std::map <std::vector<sample_id_t>, size_t> out;

	for (const auto& p : patterns) {
		++out[p];
	}
	return out;
}


/****************************************************************************************************************************************************/


FastKmerDb::FastKmerDb() : kmers2patternIds((unsigned long long) - 1), dictionarySearchQueue(1), patternExtensionQueue(1) {

	patternBytes = 0;
	patterns.reserve(1024);
	patterns.push_back(pattern_t());
	num_threads = std::thread::hardware_concurrency();

	threadPatterns.resize(num_threads);
	for (auto & tp : threadPatterns) {
		tp.reserve(1024);
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
							kmers2patternIds.prefetch(prefetch_kmer);
						}
#endif

						// ***** Sprawdzanie w slowniku czy taki k-mer juz istnieje
						auto i_kmer = kmers2patternIds.find(u_kmer);
						pattern_id_t p_id;				// tu bedzie id wzorca (pattern), ktory aktualnie jest przypisany do tego k-mera

						if (i_kmer == nullptr)
						{
							//i_kmer = kmers2patternIds.insert(u_kmer, 0);
							//p_id = 0;						// Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0

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

					(*task.num_existing_kmers)[task.block_id] = existing_id;
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

						// zliczamy odczyty z aktualnego pliku (osobnika) o tym samym wzorcu
						for (j = i + 1; j < hi; ++j) {
							if (p_id != samplePatterns[j].first) {
								break;
							}
						}
						size_t pid_count = j - i;

						if (patterns[p_id].get_num_kmers() == pid_count && !patterns[p_id].get_is_parrent()) {
							// Wzorzec mozna po prostu rozszerzyæ, bo wszystkie wskazniki do niego beda rozszerzane (realokacja tablicy id próbek we wzorcu) 
							mem -= patterns[p_id].get_bytes();
							patterns[p_id].expand(task.sample_id);
							mem += patterns[p_id].get_bytes();
						}
						else
						{
							// Trzeba wygenerowaæ nowy wzorzec (podczepiony pod wzorzec macierzysty)
							//auto *pat = new subpattern_t<sample_id_t>(*(patterns[p_id].last_subpattern), sampleId);
							//mem += pat->getMem();

							//	patterns.push_back(pattern_t{ (uint32_t)pid_count, pat });
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


FastKmerDb::~FastKmerDb() {
	dictionarySearchQueue.MarkCompleted();
	for (auto& t : dictionarySearchWorkers) {
		t.join();
	}

	patternExtensionQueue.MarkCompleted();
	for (auto& t : patternExtensionWorkers) {
		t.join();
	}
/*
	std::vector<std::thread> workers(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, tid]() {
			size_t n_patterns = patterns.size();
			size_t block_size = n_patterns / num_threads;
			size_t lo = tid * block_size;
			size_t hi = (tid == num_threads - 1) ? n_patterns : lo + block_size;

			for (size_t i = lo; i < hi; ++i) {
				patterns[i].release();
			}
		});
	}

	for (auto& t : workers) {
		t.join();
	}
	*/
}





// Przetwarza pojedyncza baze KMC
sample_id_t FastKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers)
{
#ifdef USE_PREFETCH
	const size_t prefetch_dist = 48;
#endif
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers);

	size_t n_kmers = kmers.size();
	
	LOG_DEBUG << "Restructurizing hashtable (serial)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();

	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	kmers2patternIds.reserve_for_additional(n_kmers);
	
	samplePatterns.resize(n_kmers);
#ifdef USE_RADULS
	tmp_samplePatterns.resize(n_kmers);
#endif

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
		dictionarySearchQueue.Push(DictionarySearchTask{ tid, num_threads, &kmers, &num_existing_kmers });
	}
	// wait for the task to complete
	semaphore.waitForZero();
	hashtableFindTime += std::chrono::high_resolution_clock::now() - start;
	
	// add kmers to hashtable sequentially
	LOG_DEBUG << "Adding kmers (serial)..." << endl;
	start = std::chrono::high_resolution_clock::now();

/*	for (int tid = 0; tid < num_threads; ++tid) {
		size_t n_kmers = kmers.size();
		size_t block = n_kmers / num_threads;
		size_t lo = tid * block;
		size_t hi = (tid == num_threads - 1) ? n_kmers : lo + block;

		for (size_t i = num_existing_kmers[tid]; i < hi; ++i) {
			auto i_kmer = kmers2patternIds.insert(samplePatterns[i].first, 0); // Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0				
			samplePatterns[i].first = 0;
			samplePatterns[i].second = i_kmer;
		}
	}*/

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
		auto i_kmer = kmers2patternIds.insert(samplePatterns[i].first, 0); // Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0
		samplePatterns[i].first = 0;
		samplePatterns[i].second = i_kmer;
	}

	hashtableAddTime += std::chrono::high_resolution_clock::now() - start;
	
	LOG_DEBUG << "Sorting (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();

#ifdef USE_RADULS
	uint32_t raduls_key_size;
	size_t raduls_tmp = patterns.size();
	for (raduls_key_size = 1; raduls_tmp >= 256; raduls_tmp >>= 8)
		++raduls_key_size;

	ParallelSort(samplePatterns.data(), samplePatterns.size(), tmp_samplePatterns.data(), sizeof(std::pair<pattern_id_t, pattern_id_t*>), raduls_key_size, std::thread::hardware_concurrency());
//	raduls::RadixSortMSD(reinterpret_cast<uint8_t*>(samplePatterns.data()), reinterpret_cast<uint8_t*>(tmp_samplePatterns.data()),
//		samplePatterns.size(), sizeof(std::pair<pattern_id_t, pattern_id_t*>), raduls_key_size, std::thread::hardware_concurrency());
	if (raduls_key_size & 1)
		samplePatterns.swap(tmp_samplePatterns); 
#else
	ParallelSort(samplePatterns.data(), samplePatterns.size(), nullptr, 0, 0, std::thread::hardware_concurrency());
#endif
	
	sortTime += std::chrono::high_resolution_clock::now() - start;
	std::cout << "sort time: " << sortTime.count() << "  " << samplePatterns.size() << endl;


	LOG_DEBUG << "Extending kmers (parallel)..." << endl;
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

	LOG_DEBUG << "Inserting kmers (serial)..." << endl;
	for (int tid = 0; tid < threads.size(); ++tid) {
		patternBytes += threadBytes[tid];
	}

	extensionTime += std::chrono::high_resolution_clock::now() - start;
	
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


void FastKmerDb::mapKmers2Samples(uint64_t kmer, std::vector<sample_id_t>& samples) const {
	
	// find corresponding pattern id
	auto p_id = kmers2patternIds.find(kmer);
	const auto& p = patterns[*p_id];
	samples.resize(p.get_num_samples());
	
	p.decodeSamples(samples.data());

}



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
	file.write(reinterpret_cast<char*>(&temp), sizeof(temp));
	
	// write ht elements in portions
	size_t bufpos = 0;
	for (auto it = kmers2patternIds.cbegin(); it < kmers2patternIds.cend(); ++it) {
		if (kmers2patternIds.is_free(*it)) {
			continue;
		}

		hastableBuffer[bufpos++] = *it;
		if (bufpos == numHastableElements) {
			file.write(reinterpret_cast<char*>(&bufpos), sizeof(size_t));
			file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
			bufpos = 0;
		}
	}
	// write remaining ht elements
	file.write(reinterpret_cast<char*>(&bufpos), sizeof(size_t));
	file.write(buffer, bufpos * sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t));
	
	// write patterns in portions
	temp = patterns.size();
	file.write(reinterpret_cast<char*>(&temp), sizeof(temp));
	
	char * currentPtr = buffer;
	for (int pid = 0; pid < patterns.size(); ++pid) {
		if (currentPtr + patterns[pid].get_bytes() > buffer + ioBufferBytes) {
			size_t blockSize = currentPtr - buffer;
			file.write(reinterpret_cast<char*>(&blockSize), sizeof(size_t)); // write size of block to facilitate deserialization
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
	file.write(reinterpret_cast<char*>(&blockSize), sizeof(size_t)); // write size of block to facilitate deserialization
	file.write(buffer, blockSize);
}



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
	
	return true;
}



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



void FastKmerDb::calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix) //const 
{
	int samples_count = getSamplesCount();
	matrix.resize(samples_count);
	matrix.clear();
	
	size_t bufsize = 4000000 / sizeof(uint32_t);
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

							auto *p = rawData;

#ifdef ALL_STATS
							localAdditions += num_samples;
#endif
							__m128i _to_add = _mm_set1_epi32((int)to_add);

							int j;

							if(num_samples % 32 >= 16)
							{
								j = -16;
								goto inner_start;
							}

							for (j = 0; j + 32 <= num_samples; j += 32)
							{
								if (*p + 15 == *(p + 15))
								{
									auto _q = (__m128i*) (row + *p);

									_mm_storeu_si128(_q, _mm_add_epi32(_mm_loadu_si128(_q), _to_add));
									_mm_storeu_si128(_q + 1, _mm_add_epi32(_mm_loadu_si128(_q + 1), _to_add));
									_mm_storeu_si128(_q + 2, _mm_add_epi32(_mm_loadu_si128(_q + 2), _to_add));
									_mm_storeu_si128(_q + 3, _mm_add_epi32(_mm_loadu_si128(_q + 3), _to_add));

									p += 16;
								}
								else
								{
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
								}

inner_start:
								if (*p + 15 == *(p + 15))
								{
									auto _q = (__m128i*) (row + *p);

									_mm_storeu_si128(_q, _mm_add_epi32(_mm_loadu_si128(_q), _to_add));
									_mm_storeu_si128(_q + 1, _mm_add_epi32(_mm_loadu_si128(_q + 1), _to_add));
									_mm_storeu_si128(_q + 2, _mm_add_epi32(_mm_loadu_si128(_q + 2), _to_add));
									_mm_storeu_si128(_q + 3, _mm_add_epi32(_mm_loadu_si128(_q + 3), _to_add));

									p += 16;
								}
								else
								{
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
									row[*p++] += to_add;
								}
							}
							num_samples -= j;

							switch (num_samples % 16)
							{
							case 15:	row[*p++] += to_add;
							case 14:	row[*p++] += to_add;
							case 13:	row[*p++] += to_add;
							case 12:	row[*p++] += to_add;
							case 11:	row[*p++] += to_add;
							case 10:	row[*p++] += to_add;
							case 9:		row[*p++] += to_add;
							case 8:		row[*p++] += to_add;
							case 7:		row[*p++] += to_add;
							case 6:		row[*p++] += to_add;
							case 5:		row[*p++] += to_add;
							case 4:		row[*p++] += to_add;
							case 3:		row[*p++] += to_add;
							case 2:		row[*p++] += to_add;
							case 1:		row[*p++] += to_add;
							}

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
//		if (pid > 2000000)
//			break;
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
		auto &sum_hist_sample_ids = hist_sample_ids[num_threads];

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

#if 0
void FastKmerDb::calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix) const {
	matrix.resize(getSamplesCount());
	matrix.clear();

	size_t bufsize = 8000000 / sizeof(uint32_t);
	std::vector<uint32_t> patternsBuffer(bufsize);
	uint32_t* currentPtr;

	std::vector<uint32_t*> rawPatterns(bufsize);
	std::vector<std::pair<sample_id_t, uint32_t>> sample2pattern(bufsize);

	std::vector<std::thread> workers(num_threads);

	// process all patterns in blocks determined by buffer size
	cout << endl;
	for (int pid = 0; pid < patterns.size(); ) {
		/*		if (pid > 2000000)
		break;*/
		cout << pid << " of " << patterns.size() << "\r";
		fflush(stdout);

		int first_pid = pid;
		currentPtr = patternsBuffer.data();
		size_t samplesCount = 0;

		// unpack as long as there is enough memory
		while (pid < patterns.size() && currentPtr + patterns[pid].get_num_samples() < patternsBuffer.data() + bufsize) {
			const auto& pattern = patterns[pid];

			if (pattern.get_num_kmers() > 0) {

				currentPtr += pattern.get_num_samples();
				uint32_t* out = currentPtr;		// start from the end
				samplesCount += pattern.get_num_samples();

				// decode all samples from pattern and its parents
				int64_t current_id = pid;
				while (current_id >= 0) {
					const auto& cur = patterns[current_id];
					out -= cur.get_num_local_samples();
					cur.decodeSamples(out);

					current_id = cur.get_parent_id();
				}
				rawPatterns[pid - first_pid] = out; // begin of unpacked pattern
			}
			else
				rawPatterns[pid - first_pid] = nullptr;

			++pid;
		}

		int last_pid = pid;

		// generate sample to pattern mapping
		sample2pattern.resize(samplesCount);
		int pair_id = 0;
		for (int pid = first_pid; pid < last_pid; ++pid) {
			const auto& pattern = patterns[pid];
			uint32_t* rawData = rawPatterns[pid - first_pid];

			if (pattern.get_num_kmers() > 0)
			{
				int num_samples = pattern.get_num_samples();
				for (int j = 0; j < num_samples; ++j) {
					sample2pattern[pair_id].first = rawData[j];
					sample2pattern[pair_id++].second = pid - first_pid;
				}
			}
		}

		// sort mapping wrt both elements
#ifdef WIN32
		concurrency::parallel_sort(sample2pattern.begin(), sample2pattern.end());
#else
		__gnu_parallel::sort(sample2pattern.begin(), sample2pattern.end());
#endif

		// determine ranges of blocks processed by threads 
		int no_ranges = num_threads * 16;
		std::vector<size_t> workerRanges(no_ranges + 1, sample2pattern.size());
		workerRanges[0] = 0;
		size_t workerBlock = sample2pattern.size() / no_ranges;

		auto currentIndex = workerBlock;

		for (int rid = 0; rid < no_ranges - 1; ++rid) {
			// make sure no single row of matrix is updated by multiple workers
			auto it = std::upper_bound(
				sample2pattern.begin() + currentIndex,
				sample2pattern.end(),
				*(sample2pattern.begin() + currentIndex - 1),
				[](auto x, auto y) {return x.first < y.first; });	// Necessary, because we are looking for end of sample data

			size_t range = it - sample2pattern.begin();

			workerRanges[rid + 1] = range;
			currentIndex = range + workerBlock;

			if (currentIndex >= sample2pattern.size()) {
				break;
			}
		}

		// this should never happen
		if (workerRanges[no_ranges] != sample2pattern.size()) {
			throw std::runtime_error("ERROR in FastKmerDb::calculateSimilarity(): Invalid ranges");
		}

		CRegisteringQueue<int> tasks_queue(1);
		for (int i = no_ranges - 1; i >= 0; --i)
			tasks_queue.Push(i);
		tasks_queue.MarkCompleted();

		// increment array elements in threads
		for (int tid = 0; tid < num_threads; ++tid) {
			workers[tid] = std::thread([&sample2pattern, &rawPatterns, &workerRanges, &matrix, this, &tasks_queue, first_pid] {
				// each worker processes its own block
				while (!tasks_queue.IsCompleted())
				{
					int range_id;
					if (tasks_queue.Pop(range_id))
					{
						for (int id = workerRanges[range_id]; id < workerRanges[range_id + 1]; ) {
							int Si = sample2pattern[id].first;
							uint32_t *row = matrix[Si];
							while (id < workerRanges[range_id + 1] && sample2pattern[id].first == Si) {
								int local_pid = sample2pattern[id].second;
								const auto& pattern = patterns[local_pid + first_pid];

								uint32_t* rawData = rawPatterns[local_pid];
								int num_samples = pattern.get_num_samples();
								uint32_t to_add = pattern.get_num_kmers();
								/*								for (int j = 0; j < num_samples; ++j) {
								uint32_t Sj = rawData[j];
								if (Sj < Si)
								row[Sj] += to_add;
								else
								break;
								}*/

								if (num_samples < 15)
								{
									int j = 0;
									uint32_t Sj;
									switch (num_samples % 4)
									{
									case 3:
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
									case 2:
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
									case 1:
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
									}
									for (; j < num_samples && rawData[j] < Si;)
									{
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
										Sj = rawData[j++];	if (Sj < Si)	row[Sj] += to_add;
									}
								}
								else
								{
									num_samples = lower_bound(rawData, rawData + num_samples, Si) - rawData;
									auto *p = rawData;

									switch (num_samples % 8)
									{
									case 7:	row[*p++] += to_add;
									case 6:	row[*p++] += to_add;
									case 5:	row[*p++] += to_add;
									case 4:	row[*p++] += to_add;
									case 3:	row[*p++] += to_add;
									case 2:	row[*p++] += to_add;
									case 1:	row[*p++] += to_add;
									}
									for (int j = num_samples % 8; j < num_samples; j += 8)
									{
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
										row[*p++] += to_add;
									}
								}

								++id;
							}
						}
					}
				}
			});
		}

		for (auto & w : workers) {
			w.join();
		}
	}
}
#endif

/*
#define BUFFERED_ARRAY

void FastKmerDb::calculateSimilarity(LowerTriangularMatrix<uint32_t>& matrix) const {
	if(getSamplesCount() < 12000)
		calculateSimilarityDirect(matrix);
	else
#ifdef BUFFERED_ARRAY
	calculateSimilarityBuffered(matrix);
#else
	calculateSimilarityDirect(matrix);
#endif
}


*/
void FastKmerDb::calculateSimilarityDirect(LowerTriangularMatrix<uint32_t>& matrix) const {
	matrix.resize(getSamplesCount());
	matrix.clear();
	
	std::vector<std::thread> workers(num_threads);
	std::atomic<uint64_t> numAdditions(0);

	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, tid, &matrix, &numAdditions]() {
			uint64_t localAdditions = 0;
			std::vector<uint32_t> rawData(getSamplesCount());

			for (int pid = 1; pid < patterns.size(); ++pid) {
				const auto& pattern = patterns[pid];
				uint32_t* out = rawData.data() + pattern.get_num_samples(); // start from the end
														// decode all samples

				if (pattern.get_num_kmers() == 0)
					continue;

				int64_t current_id = pid;
				while (current_id >= 0) {
					const auto& cur = patterns[current_id];
					out -= cur.get_num_local_samples();
					cur.decodeSamples(out);

					current_id = cur.get_parent_id();
				}

				for (int i = 0; i < pattern.get_num_samples(); ++i) {
					uint32_t Si = rawData[i];
					uint32_t key = Si % (2 * num_threads);

					if (key == tid || key == (2 * num_threads - 1 - tid)) {
						uint32_t * row = matrix[Si];
						for (int j = 0; j < i; ++j) {
							uint32_t Sj = rawData[j];
							row[Sj] += pattern.get_num_kmers();
#ifdef ALL_STATS
							++localAdditions;
#endif
						}
					}
				}
			}
			numAdditions.fetch_add(localAdditions);
		});
	}
	
	for (auto & w : workers) {
		w.join();
	}

#ifdef ALL_STATS
	cout << "Number of additions:" << numAdditions << endl;
#endif
	
}

void FastKmerDb::calculateSimilarityBuffered(LowerTriangularMatrix<uint32_t>& matrix) const {
	matrix.resize(getSamplesCount());
	matrix.clear();

	std::vector<std::thread> workers(num_threads);
	std::atomic<uint64_t> numAdditions(0);

	vector<ArrayBuffer> mat_buf(getSamplesCount());
	for (size_t i = 0; i < getSamplesCount(); ++i)
		mat_buf[i].Assign(matrix[i], 1 << 15);

	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, tid, &matrix, &numAdditions, &mat_buf]() {
			uint64_t localAdditions = 0;
			std::vector<uint32_t> rawData(getSamplesCount());
			std::vector<uint32_t> v_samples;

			for (int pid = 1; pid < patterns.size(); ++pid) {
				const auto& pattern = patterns[pid];
				if (pattern.get_num_kmers() == 0)
					continue;

				uint32_t* out = rawData.data() + pattern.get_num_samples(); // start from the end
																			// decode all samples

				int64_t current_id = pid;
				while (current_id >= 0) {
					const auto& cur = patterns[current_id];
					out -= cur.get_num_local_samples();
					cur.decodeSamples(out);

					current_id = cur.get_parent_id();
				}

				v_samples.clear();
				for (int i = 0; i < pattern.get_num_samples(); ++i) {
					uint32_t Si = rawData[i];
					uint32_t key = Si % (2 * num_threads);

					if (key == tid || key == (2 * num_threads - 1 - tid))
						v_samples.push_back(i);
				}

				uint32_t n_samples = v_samples.size();
//				for (int i = 0; i < pattern.get_num_samples() - 1; ++i) {
				for(uint32_t ii  = 0; ii < n_samples; ++ii)
				{
					if (ii + 2 < n_samples)
					{
						_mm_prefetch((const char*)(rawData.data() + v_samples[ii + 1]), _MM_HINT_T0);
						mat_buf[rawData[v_samples[ii + 2]]].Prefetch();
					}

					uint32_t i = v_samples[ii];
					uint32_t Si = rawData[i];

					ArrayBuffer &row_buf = mat_buf[Si];
					row_buf.SetCounter(pattern.get_num_kmers());

/*					for (int j = 0; j < i; ++j) {
						row_buf.Push(rawData[j]);*/
					auto *end_p = rawData.data() + i;
					for (auto p = rawData.data(); p != end_p; ++p){
						row_buf.Push(*p);
#ifdef ALL_STATS
						++localAdditions;
#endif
					}
				}
			}
			numAdditions.fetch_add(localAdditions);

			for (size_t i = 0; i < getSamplesCount(); ++i)
			{
				uint32_t key = i % (2 * num_threads);
				if (key == tid || key == (2 * num_threads - 1 - tid))
					mat_buf[i].Finish();
			}
		});
	}

	for (auto & w : workers) {
		w.join();
	}

#ifdef ALL_STATS
	cout << "Number of additions:" << numAdditions << endl;
#endif

}

void  FastKmerDb::calculateSimilarity(const FastKmerDb& sampleDb, std::vector<uint32_t>& similarities) const {
	similarities.resize(this->getSamplesCount(), 0);

	std::vector<std::vector<uint32_t>> localSimilarities(num_threads, std::vector<uint32_t>(this->getSamplesCount()));

	std::unordered_map<pattern_id_t, int32_t> patterns2count;

	std::chrono::duration<double> dt;
	auto start = std::chrono::high_resolution_clock::now();

	// iterate over kmers in analyzed sample
	for (auto it = sampleDb.kmers2patternIds.cbegin(); it < sampleDb.kmers2patternIds.cend(); ++it) {
		if (sampleDb.kmers2patternIds.is_free(*it)) {
			continue;
		}

		// check if kmer exists in a database
		auto entry = kmers2patternIds.find(it->key);

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
	for (const auto& entry : patterns2count) {
		patterns2countVector[i++] = entry;
	}
	patterns2count.clear();

	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "Pattern listing time: " << dt.count() << endl;

	start = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> workers(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, tid, &patterns2countVector, &localSimilarities]() {
			std::vector<uint32_t> samples(this->getSamplesCount());
			
			size_t n_patterns = patterns2countVector.size();
			size_t block_size = n_patterns / num_threads;
			size_t lo = tid * block_size;
			size_t hi = (tid == num_threads - 1) ? n_patterns : lo + block_size;
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
	cout << "Pattern unpacking time: " << dt.count() << endl;
}
