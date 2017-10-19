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

#ifndef WIN32
#include <parallel/algorithm>
#endif

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"

using namespace std;


bool AbstractKmerDb::extractKmers(const string &filename, std::vector<uint64_t>& kmers) {
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



FastKmerDb::FastKmerDb() {
	hm_kmer_dict.set_special_keys((unsigned long long) - 1, (unsigned long long) - 2);
	v_kmer_patterns.push_back(pattern_t{ 0, new vid_t() });
}


// Przetwarza pojedyncza baze KMC
void FastKmerDb::addKmers(sample_id_t sampleId, const std::vector<uint64_t>& kmers)
{
	uint64_t u_kmer;
	uint64_t prefetch_kmer;
	const size_t prefetch_dist = 48;
	auto n_kmers = kmers.size();
	
	v_current_file_pids.clear();
	
	// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	hm_kmer_dict.reserve_for_additional(n_kmers);

	for (size_t i = 0; i < n_kmers; ++i)
	{
		u_kmer = kmers[i];
		if (i + prefetch_dist < n_kmers)
		{
			prefetch_kmer = kmers[i + prefetch_dist];
			hm_kmer_dict.prefetch(prefetch_kmer);
		}

		// ***** Sprawdzanie w slowniku czy taki k-mer juz istnieje
		auto i_kmer = hm_kmer_dict.find(u_kmer);
		uint64_t p_id;				// tu bedzie id wzorca (pattern), ktory aktualnie jest przypisany do tego k-mera

		if (i_kmer == nullptr)
		{
			i_kmer = hm_kmer_dict.insert(u_kmer, 0);
			p_id = 0;						// Pierwsze wyst. k-mera, wiec przypisujemy taki sztuczny wzorzec 0
		}
		else
			p_id = *i_kmer;

		v_current_file_pids.push_back(make_pair(p_id, i_kmer));
//		hm_current_file_pids.insert(make_pair(p_id, i_kmer));
	}

#ifdef WIN32
	sort(v_current_file_pids.begin(), v_current_file_pids.end());
#else
	__gnu_parallel::sort(v_current_file_pids.begin(), v_current_file_pids.end());
#endif

	for (size_t i = 0; i < n_kmers;)
	{
		size_t j;
		auto p_id = v_current_file_pids[i].first;

		// zliczamy odczyty z aktualnego pliku (osobnika) o tym samym wzorcu
		for (j = i + 1; j < n_kmers; ++j)
			if (p_id != v_current_file_pids[j].first)
				break;
		size_t pid_count = j - i; 

		if (v_kmer_patterns[p_id].num_occ == pid_count && !v_kmer_patterns[p_id].vid->get_is_parrent())
		{
			// Wzorzec mozna po prostu rozszerzyæ, bo wszystkie wskazniki do niego beda rozszerzane (realokacja tablicy id próbek we wzorcu) 
			v_kmer_patterns[p_id].vid->expand(sampleId);
		}
		else
		{
			// Trzeba wygenerowaæ nowy wzorzec (podczepiony pod wzorzec macierzysty)
			vid_t *pat = new vid_t(*(v_kmer_patterns[p_id].vid), sampleId);

			v_kmer_patterns.push_back(pattern_t{ (uint32_t) pid_count, pat });
			if(p_id)
				v_kmer_patterns[p_id].num_occ -= pid_count;

			uint32_t new_p_id = v_kmer_patterns.size() - 1;

			for (size_t k = i; k < j; ++k)
				*(v_current_file_pids[k].second) = new_p_id;
		}

		i = j;
	}
}


void FastKmerDb::kmer2samples(uint64_t kmer, std::vector<sample_id_t>& samples) {
	
	// find corresponding pattern id
	auto p_id = hm_kmer_dict.find(kmer);
	samples.clear();
	if (p_id != nullptr) { 
		auto subpattern = v_kmer_patterns[*p_id].vid;
		
		do {
			
		} while (subpattern != nullptr);
	}


}