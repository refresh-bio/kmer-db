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


#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "log.h"

#ifdef USE_RADULS
#include "raduls/raduls.h"
#endif

using namespace std;

const size_t FastKmerDb::ioBufferBytes = (2 << 29); //512MB buffer 
//const size_t FastKmerDb::ioBufferBytes = 16000000; //16MB buffer 

/****************************************************************************************************************************************************/

bool AbstractKmerDb::loadKmers(const string &filename, std::vector<kmer_t>& kmers) {
	CKMCFile kmc_file;
	uint32_t counter;

	if (!kmc_file.OpenForListing(filename))
		return false;

	uint32 _kmer_length;
	uint32 _mode;
	uint32 _counter_size;
	uint32 _lut_prefix_length;
	uint32 _signature_len;
	uint32 _min_count;
	uint64 _max_count;
	uint64 _total_kmers;

	kmc_file.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	CKmerAPI kmer(_kmer_length);

	uint64_t u_kmer;
	vector<uint64> tmp;

	// Wczytuje wszystkie k-mery z pliku do wektora, zeby pozniej moc robic prefetcha
	kmers.clear();
	
	while (!kmc_file.Eof())
	{
		if (!kmc_file.ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		kmers.push_back(u_kmer);
	}

	kmc_file.Close();

	return true;
}


/****************************************************************************************************************************************************/

sample_id_t NaiveKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers);
	
	uint64_t u_kmer;
	uint64_t prefetch_kmer;
	const size_t prefetch_dist = 48;
	auto n_kmers = kmers.size();
	
	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	kmers2patternIds.reserve_for_additional(n_kmers);
	patterns.reserve(patterns.size() + n_kmers);

	for (size_t i = 0; i < n_kmers; ++i) {
		u_kmer = kmers[i];
		if (i + prefetch_dist < n_kmers) {
			prefetch_kmer = kmers[i + prefetch_dist];
			kmers2patternIds.prefetch(prefetch_kmer);
		}

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
	patterns.push_back(pattern_t<sample_id_t>());
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
					kmer_t prefetch_kmer;
					const size_t prefetch_dist = 48;

					for (size_t i = lo; i < hi; ++i)
					{
						u_kmer = kmers[i];
						if (i + prefetch_dist < hi)
						{
							prefetch_kmer = kmers[i + prefetch_dist];
							kmers2patternIds.prefetch(prefetch_kmer);
						}

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

							threadPatterns[task.block_id].emplace_back(local_pid, pattern_t<sample_id_t>(patterns[p_id], p_id, task.sample_id, (uint32_t)pid_count));
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
}





// Przetwarza pojedyncza baze KMC
sample_id_t FastKmerDb::addKmers(std::string sampleName, const std::vector<kmer_t>& kmers)
{
	sample_id_t sampleId = AbstractKmerDb::addKmers(sampleName, kmers);

	size_t n_kmers = kmers.size();
	
	LOG_DEBUG << "Restructurizing hashtable (serial)..." << endl;
	auto start = std::chrono::high_resolution_clock::now();

	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	kmers2patternIds.reserve_for_additional(n_kmers);
	
	samplePatterns.resize(n_kmers);
#ifdef RADULS
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

	for (int tid = 0; tid < num_threads; ++tid) {
		size_t n_kmers = kmers.size();
		size_t block = n_kmers / num_threads;
		size_t lo = tid * block;
		size_t hi = (tid == num_threads - 1) ? n_kmers : lo + block;

		for (size_t i = num_existing_kmers[tid]; i < hi; ++i) {
			auto i_kmer = kmers2patternIds.insert(samplePatterns[i].first, 0); // Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0				
			samplePatterns[i].first = 0;
			samplePatterns[i].second = i_kmer;
		}
	}

	hashtableAddTime += std::chrono::high_resolution_clock::now() - start;
	
	LOG_DEBUG << "Sorting (parallel)..." << endl;
	start = std::chrono::high_resolution_clock::now();
	auto pid_comparer = [](const std::pair<pattern_id_t, pattern_id_t*>& a, const std::pair<pattern_id_t, pattern_id_t*>& b)->bool {
		return a.first < b.first;
	};
#ifdef USE_RADULS
	raduls::RadixSortMSD(reinterpret_cast<uint8_t*>(samplePatterns.data()), reinterpret_cast<uint8_t*>(tmp_samplePatterns.data()), 
//		samplePatterns.size(), sizeof(std::pair<pattern_id_t, pattern_id_t*>), sizeof(size_t), std::thread::hardware_concurrency());
		samplePatterns.size(), sizeof(std::pair<pattern_id_t, pattern_id_t*>), 5, 12);
#else
#ifdef WIN32
	concurrency::parallel_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#else
	__gnu_parallel::sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#endif
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

	for (int tid = 0; tid < num_threads; ++tid) {
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
	samples.clear();
	if (p_id != nullptr) { 
		const pattern_t<sample_id_t>* subpattern = &patterns[*p_id];
		samples.resize(subpattern->get_num_samples());
		int j = samples.size() - 1;

		// collect in reversed order
		int64_t cur_id = *p_id;

		do {
			subpattern = &patterns[cur_id];

			for (int i = subpattern->get_num_local_samples() - 1; i >= 0 ; --i, --j) {
				samples[j] = (*subpattern)[i];
			}
			
			cur_id = subpattern->get_parent_id(); // return -1 when no parrent
			
		} while (cur_id >= 0);
	}
}



void FastKmerDb::serialize(std::ofstream& file) const {
	
	size_t numHastableElements = ioBufferBytes / sizeof(hash_map_lp<kmer_t, pattern_id_t>::item_t);
	std::vector <hash_map_lp<kmer_t, pattern_id_t>::item_t> hastableBuffer(numHastableElements);
	char* buffer = reinterpret_cast<char*>(hastableBuffer.data());

	// store number of samples
	size_t temp = getSamplesCount();
	file.write(reinterpret_cast<const char*>(&temp), sizeof(temp));

	// store sample names
	for (const string& s : sampleNames) {
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

	// load sample names
	size_t temp;
	file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	sampleNames.resize(temp);

	for (string& s : sampleNames) {
		file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
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

void FastKmerDb::calculateSimilarity(Array<uint32_t>& matrix) const {
	matrix.resize(getSamplesCount(), getSamplesCount());
	matrix.clear();

	std::vector<uint32_t> rawData(getSamplesCount());
	
	for (int pid = 1; pid < patterns.size(); ++pid) {
		const auto& pattern = patterns[pid];
		uint32_t* pos = rawData.data();

		// decode all samples
		int64_t current_id = pid;
		while (current_id >= 0) {
			const auto& cur = patterns[current_id];
			cur.decodeSamples(pos);
			pos += cur.get_num_local_samples();

			current_id = cur.get_parent_id();
		}

		for (int i = 0; i < pattern.get_num_samples() - 1; ++i) {
			sample_id_t Si = rawData[i];

			for (int j = i + 1; j < pattern.get_num_samples(); ++j) {
				sample_id_t Sj = rawData[j];
				sample_id_t first, last;
				if (Si < Sj) {
					first = Si;
					last = Sj;
				}
				else {
					first = Sj;
					last = Si;
				}

				matrix[first][last] += pattern.get_num_kmers();
				//matrix[Sj][Si] += pattern.get_num_kmers();
			}
		}
	}
}

void  FastKmerDb::calculateSimilarity(const FastKmerDb& sampleDb, std::vector<uint32_t>& similarities) const {
		similarities.resize(this->getSamplesCount());

		// iterate over kmers in analyzed sample
		for (auto it = sampleDb.kmers2patternIds.cbegin(); it < sampleDb.kmers2patternIds.cend(); ++it) {
			if (sampleDb.kmers2patternIds.is_free(*it)) {
				continue;
			}

			// check if kmer exists in a database
			auto entry = kmers2patternIds.find(it->key);

			if (entry != nullptr) {
				auto pid = *entry;
				while (pid >= 0) {
					const auto& pat = patterns[pid];

					// iterate over elements in the subpattern
					for (int i = pat.get_num_local_samples() - 1; i >= 0; --i) {
						sample_id_t Si = pat.get_data()[i];
						++similarities[Si];
					}

					pid = pat.get_parent_id();
				}
			}
		}
	}