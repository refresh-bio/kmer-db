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

#ifndef WIN32
#include <parallel/algorithm>
#endif

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"

using namespace std;

/****************************************************************************************************************************************************/

bool AbstractKmerDb::extractKmers(const string &filename, std::vector<kmer_t>& kmers) {
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

void NaiveKmerDb::addKmers(sample_id_t sampleId, const std::vector<kmer_t>& kmers)
{
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

		mem_pattern_desc += sizeof(sample_id_t);
	}
}


void NaiveKmerDb::mapKmers2Samples(kmer_t kmer, std::vector<sample_id_t>& samples) {
	samples.clear();
	auto v = kmers2patternIds.find(kmer);

	if (v != nullptr) {
		samples = patterns[*v];
	}
}



/****************************************************************************************************************************************************/


FastKmerDb::FastKmerDb() : kmers2patternIds((unsigned long long) - 1, (unsigned long long) - 2) {

	mem_pattern_desc = 0;
	patterns.push_back(pattern_t{ 0, new subpattern_t<sample_id_t>() });
	
	/*
	samplePatterns.reserve(2 << 23); // 16M elements
	
	for (int i = 0; i < threadsSamplePatterns.size(); ++i) {
		threadsSamplePatterns[i].reserve(samplePatterns.size() / threadsSamplePatterns.size());
	}*/
	
}


// Przetwarza pojedyncza baze KMC
void FastKmerDb::addKmers(sample_id_t sampleId, const std::vector<kmer_t>& kmers)
{
	size_t n_kmers = kmers.size();
	
	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	kmers2patternIds.reserve_for_additional(n_kmers);
	std::mutex lock;
	
	samplePatterns.resize(n_kmers);
	int n_threads = std::thread::hardware_concurrency();
	std::vector<size_t> threadsExistingSamples(n_threads);

	std::vector<std::thread> threads(n_threads);
	for (int tid = 0; tid < threads.size(); ++tid) {
		threads[tid] = std::thread([this, &kmers, tid, n_threads, &threadsExistingSamples]() {
			
			size_t n_kmers = kmers.size();
			size_t block = n_kmers / n_threads;
			size_t lo = tid * block;
			size_t hi = (tid == n_threads - 1) ? n_kmers : lo + block;

		
			size_t existing_sample_id = lo;
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
				pid_t p_id;				// tu bedzie id wzorca (pattern), ktory aktualnie jest przypisany do tego k-mera

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

					samplePatterns[existing_sample_id].first = p_id;
					samplePatterns[existing_sample_id].second = i_kmer;
					existing_sample_id++;
				}
			}
			threadsExistingSamples[tid] = existing_sample_id;
			cout << "Thread " << tid << ", existing: " << existing_sample_id - lo << ", to add: " << hi - existing_sample_id << endl;
		});
	}

	size_t n_patterns = 0;
	for (int tid = 0; tid < threads.size(); ++tid) {
		threads[tid].join();
	}

	// add kmers to hashtable sequentially
	for (int tid = 0; tid < threads.size(); ++tid) {
		size_t n_kmers = kmers.size();
		size_t block = n_kmers / n_threads;
		size_t lo = tid * block;
		size_t hi = (tid == n_threads - 1) ? n_kmers : lo + block;

		for (size_t i = threadsExistingSamples[tid]; i < hi; ++i) {
			auto i_kmer = kmers2patternIds.insert(samplePatterns[i].first, 0); // Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0				
			samplePatterns[i].first = 0;
			samplePatterns[i].second = i_kmer;
		}
	}


		
#ifdef WIN32
	sort(samplePatterns.begin(), samplePatterns.end(), [](std::pair<pid_t, pid_t*>& a, std::pair<pid_t, pid_t*>& b)->bool {
		return a.first < b.first;
	});
#else
	__gnu_parallel::sort(samplePatterns.begin(), samplePatterns.end(), [](std::pair<pid_t, pid_t*>& a, std::pair<pid_t, pid_t*>& b)->bool {
		return a.first < b.first;
	});
#endif

	std::atomic<int> new_pid = patterns.size();

	for (size_t i = 0; i < n_kmers;)
	{
		size_t j;
		auto p_id = samplePatterns[i].first;

		// zliczamy odczyty z aktualnego pliku (osobnika) o tym samym wzorcu
		for (j = i + 1; j < n_kmers; ++j)
			if (p_id != samplePatterns[j].first)
				break;
		size_t pid_count = j - i; 

		if (patterns[p_id].num_occ == pid_count && !patterns[p_id].last_subpattern->get_is_parrent())
		{
			// Wzorzec mozna po prostu rozszerzyæ, bo wszystkie wskazniki do niego beda rozszerzane (realokacja tablicy id próbek we wzorcu) 
			mem_pattern_desc -= patterns[p_id].last_subpattern->getMem();
			patterns[p_id].last_subpattern->expand(sampleId);
			mem_pattern_desc += patterns[p_id].last_subpattern->getMem();
		}
		else
		{
			// Trzeba wygenerowaæ nowy wzorzec (podczepiony pod wzorzec macierzysty)
			auto *pat = new subpattern_t<sample_id_t>(*(patterns[p_id].last_subpattern), sampleId);

			mem_pattern_desc += pat->getMem();

			patterns.push_back(pattern_t{ (uint32_t) pid_count, pat });
			if (p_id) {
				patterns[p_id].num_occ -= pid_count;
			}

			for (size_t k = i; k < j; ++k) {
				*(samplePatterns[k].second) = new_pid;
			}
			++new_pid;
		}

		i = j;
	}
}


void FastKmerDb::mapKmers2Samples(uint64_t kmer, std::vector<sample_id_t>& samples) {
	
	// find corresponding pattern id
	auto p_id = kmers2patternIds.find(kmer);
	samples.clear();
	if (p_id != nullptr) { 
		const subpattern_t<sample_id_t>* subpattern = patterns[*p_id].last_subpattern;
		samples.resize(subpattern->get_num_samples());
		int j = samples.size() - 1;

		// collect in reversed order
		do {
			for (int i = subpattern->get_num_local_samples() - 1; i >= 0 ; --i, --j) {
				samples[j] = (*subpattern)[i];
			}
			subpattern = subpattern->get_parent();

		} while (subpattern != nullptr);
	}


}