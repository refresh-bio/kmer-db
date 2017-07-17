// kmer_db.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

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

#include "cvector.h"
#include "hashmap_lp.h"
#include "kmc_api/kmc_file.h"

using namespace std;

//#define PLAIN_VECTORS

vector<string> kmc_file_list;

// K-mer database structures
hash_map<uint64_t, uint32_t> hm_kmer_dict;
//hash_map<uint64_t, uint64_t> hm_kmer_dict;


#ifdef PLAIN_VECTORS
//typedef const_vector<uint16_t> vid_t;
typedef fake_const_vector<uint16_t> vid_t;

typedef struct
{
	uint32_t num_occ;
	vid_t vid;
	uint32_t mapped_to_id;
} pattern_t;
#else
typedef pattern_desc<uint16_t> vid_t;

typedef struct
{
	uint32_t num_occ;
	vid_t *vid;
//	uint32_t mapped_to_id;
} pattern_t;

vector<pair<uint64_t, uint32_t*>> v_current_file_pids;

unordered_multimap<uint32_t, uint32_t*> hm_current_file_pids;
#endif


//vector<pair<uint32_t, vid_t>> v_kmer_patterns;					// (liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
vector<pattern_t> v_kmer_patterns;					// (liczba wystapien wzorca w kmer_dict, wektor id osobnikow zawierajacych k - mer)
stack<uint32_t, vector<uint32_t>> free_pattern_ids;

// Id wzorcow, ktore dla biezacego genomu zostaly juz uzyte
stack<uint32_t, vector<uint32_t>> s_used_mappings;

vector<uint64_t> v_current_file_kmers;
vector<uint32_t> v_decrease_counts;

uint64_t last_pattern_id = 0;
std::chrono::duration<double> file_time;

//map<uint32_t, uint32_t> pat_sizes;					// (rozmiar wzorca, liczba wyst wzorcow z takim rozmiarem)
vector<uint32_t> pat_sizes;						// rozmiar wzorca -> liczba wyst wzorcow z takim rozmiarem


// Functions
bool load_file_list(string file_list);
bool process_single_file(uint16_t current_file_id, string file_name);
void show_progress();
void store_pat_sizes(int num);

// Wczytuje plik zawierajacy liste baz KMC do przetworzenia
bool load_file_list(string file_list)
{
	char s[1024];

	kmc_file_list.clear();

	FILE *f = fopen(file_list.c_str(), "rt");
	if (!f)
		return false;

	while(!feof(f))
	{
		fscanf(f, "%s", s);
		kmc_file_list.push_back(string(s));
	}

	fclose(f);

	// Dla uproszczenia taki filtr
	// Chodzi o to, ze plik z wykazem baz KMC moze miec zdublowane nazwy plikow i nalezy wyciagnac unikalne nazwy
	sort(kmc_file_list.begin(), kmc_file_list.end());
	kmc_file_list.erase(unique(kmc_file_list.begin(), kmc_file_list.end()), kmc_file_list.end());

	return !file_list.empty();
}

#ifdef PLAIN_VECTORS
// Przetwarza pojedyncza baze KMC
bool process_single_file(uint16_t current_file_id, string file_name)
{
	CKMCFile kmc_file;
	uint32_t counter;

	if (!kmc_file.OpenForListing(file_name))
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
	uint64_t prefetch_kmer;
	const size_t prefetch_dist = 48;
	vector<uint64> tmp;

	uint64_t i = 0;

//	hm_current_mappings.clear();

//	v_decrease_counts.clear();

	// Wczytuje wszystkie k-mery z pliku do wektora, zeby pozniej moc robic prefetcha
	v_current_file_kmers.clear();

	while (!kmc_file.Eof())
	{
		if (!kmc_file.ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		v_current_file_kmers.push_back(u_kmer);
	}

/*	while (!kmc_file.Eof())
	{
		if (!kmc_file.ReadNextKmer(kmer, counter))
			break;
#ifdef WIN32
//		++i;
		if (++i % (1 << 18) == 0)
			cout << i << " of " << _total_kmers << "\r";
#endif

		// ***** Konwersja k-mera na uint64_t
		kmer.to_long(tmp);
		u_kmer = tmp.front();*/

	auto n_kmers = v_current_file_kmers.size();
	for (size_t i = 0; i < n_kmers; ++i)
	{
		u_kmer = v_current_file_kmers[i];
		if (i + prefetch_dist < n_kmers)
		{
			prefetch_kmer = v_current_file_kmers[i + prefetch_dist];
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

		// ***** Rozszerzanie wzorca zawierajacego id osobnikow 
		
		// Sprawdzanie czy juz taki wzorzec byl rozszerzany
//		auto i_mapping = hm_current_mappings.find(p_id);
//		auto i_mapping = v_current_mappings[p_id];
		auto i_mapping = v_kmer_patterns[p_id].mapped_to_id;

//		if (i_mapping == nullptr)
		if(i_mapping == 0)
		{
			// Wzorzec jeszcze nie byl rozszerzany, wiec trzeba zrobic nowy
//			vid_t pat(v_kmer_patterns[p_id].second, current_file_id);
			vid_t pat(v_kmer_patterns[p_id].vid, current_file_id);

			if (free_pattern_ids.empty())
			{
				size_t vs = v_kmer_patterns.size();
//				v_kmer_patterns.resize(vs * 2, make_pair(0, vid_t()));
				v_kmer_patterns.resize(vs * 2, pattern_t{ 0, vid_t(), 0 });
				for (uint32_t i = vs; i < vs * 2; ++i)
					free_pattern_ids.push(i);

//				v_current_mappings.resize(vs * 2, 0);
			}

			last_pattern_id = free_pattern_ids.top();
			free_pattern_ids.pop();
			
			size_t pat_size = pat.get_size();

//			v_kmer_patterns[last_pattern_id] = make_pair(0, pat);
			v_kmer_patterns[last_pattern_id] = pattern_t{ 0, pat, 0 };

//			i_mapping = hm_current_mappings.insert(p_id, last_pattern_id);
//			i_mapping = v_current_mappings[p_id] = last_pattern_id;
			i_mapping = v_kmer_patterns[p_id].mapped_to_id = last_pattern_id;

			s_used_mappings.push(p_id);

			pat_sizes[pat_size] += 1;
		}

//		*i_kmer = *i_mapping;
		*i_kmer = i_mapping;

//		++v_kmer_patterns[*i_mapping].first;
//		++v_kmer_patterns[i_mapping].first;
		++v_kmer_patterns[i_mapping].num_occ;

//		if (--(v_kmer_patterns[p_id].first) == 0)		// decrease no. of occ. of old pattern
		if (--(v_kmer_patterns[p_id].num_occ) == 0)		// decrease no. of occ. of old pattern
		{
//			pat_sizes[v_kmer_patterns[p_id].second.get_size()] -= 1;
			pat_sizes[v_kmer_patterns[p_id].vid.get_size()] -= 1;
//			v_kmer_patterns[p_id] = make_pair(0, vid_t());
			v_kmer_patterns[p_id] = pattern_t{ 0, vid_t(), 0 };
			free_pattern_ids.push(p_id);
		}
	}

	// Usuwanie mapowañ z wektora wg danych ze stosu
	while (!s_used_mappings.empty())
	{
		auto x = s_used_mappings.top();
		s_used_mappings.pop();
//		v_current_mappings[x] = 0;
		v_kmer_patterns[x].mapped_to_id = 0;
	}

	kmc_file.Close();

	return true;
}
#else
// Przetwarza pojedyncza baze KMC
bool process_single_file(uint16_t current_file_id, string file_name)
{
	CKMCFile kmc_file;
	uint32_t counter;

	if (!kmc_file.OpenForListing(file_name))
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
	uint64_t prefetch_kmer;
	const size_t prefetch_dist = 48;
	vector<uint64> tmp;

	uint64_t i = 0;

	// Wczytuje wszystkie k-mery z pliku do wektora, zeby pozniej moc robic prefetcha
	v_current_file_kmers.clear();
	v_current_file_pids.clear();
	hm_current_file_pids.clear();

	while (!kmc_file.Eof())
	{
		if (!kmc_file.ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		v_current_file_kmers.push_back(u_kmer);
	}

	auto n_kmers = v_current_file_kmers.size();
	
		// Musimy miec poprawne wskazniki do HT po wstawieniu do n_kmers elementow, a wiec nie moze sie w tym czasie zrobic restrukturyzacja
	// Jesli jest ryzyko, to niech zrobi sie wczesniej, przed wstawianiem elementow
	hm_kmer_dict.reserve_for_additional(n_kmers);

	for (size_t i = 0; i < n_kmers; ++i)
	{
		u_kmer = v_current_file_kmers[i];
		if (i + prefetch_dist < n_kmers)
		{
			prefetch_kmer = v_current_file_kmers[i + prefetch_dist];
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

		for (j = i + 1; j < n_kmers; ++j)
			if (p_id != v_current_file_pids[j].first)
				break;
		size_t pid_count = j - i;

		if (v_kmer_patterns[p_id].num_occ == pid_count && !v_kmer_patterns[p_id].vid->get_is_parrent())
		{
			// Wzorzec mozna po prostu rozszerzyæ, bo wszystkie wskazniki do niego beda rozszerzane
			v_kmer_patterns[p_id].vid->expand(current_file_id);
		}
		else
		{
			// Trzeba wygenerowaæ nowy wzorzec
			vid_t *pat = new vid_t(*(v_kmer_patterns[p_id].vid), current_file_id);

			v_kmer_patterns.push_back(pattern_t{ (uint32_t) pid_count, pat });
			if(p_id)
				v_kmer_patterns[p_id].num_occ -= pid_count;

			uint32_t new_p_id = v_kmer_patterns.size() - 1;

			for (size_t k = i; k < j; ++k)
				*(v_current_file_pids[k].second) = new_p_id;
		}

		i = j;
	}
/*	for (auto p = hm_current_file_pids.begin(); p != hm_current_file_pids.end();)
	{ 
		auto q = p;
		int pid_count = 0;
		for (; q != hm_current_file_pids.end() && p->first == q->first; ++q, ++pid_count)
			;

		auto p_id = p->first;

		if (v_kmer_patterns[p_id].num_occ == pid_count && !v_kmer_patterns[p_id].vid->get_is_parrent())
		{
			// Wzorzec mozna po prostu rozszerzyæ, bo wszystkie wskazniki do niego beda rozszerzane
			v_kmer_patterns[p_id].vid->expand(current_file_id);
		}
		else
		{
			// Trzeba wygenerowaæ nowy wzorzec
			vid_t *pat = new vid_t(*(v_kmer_patterns[p_id].vid), current_file_id);

			v_kmer_patterns.push_back(pattern_t{ (uint32_t)pid_count, pat });
			if (p_id)
				v_kmer_patterns[p_id].num_occ -= pid_count;

			uint32_t new_p_id = v_kmer_patterns.size() - 1;

			for (auto r = p; r != q; ++r)
				*(r->second) = new_p_id;
		}

		p = q;
	}*/

	kmc_file.Close();

	return true;
}
#endif
// Pokazuje stan
void show_progress()
{
	size_t tot_pat_size = 0;
	size_t num_calc = 0;			// Liczba operacji przy wyznaczaniu macierzy podobieñstwa

#ifdef PLAIN_VECTORS

	for (int i = 0; i < pat_sizes.size(); ++i)
	{
		auto x = pat_sizes[i];

		if (!x)
			continue;

		tot_pat_size += sizeof(uint32_t) * x;		// id wzorca
		tot_pat_size += sizeof(uint32_t) * x;		// liczba wyst. wzorca
		tot_pat_size += sizeof(uint16_t) * i * x;		// id elementow wzorca
		
		num_calc += i * (i - 1) / 2 * x;
	}

	cout << "*** dict: " << hm_kmer_dict.get_size()
		<< " (" << hm_kmer_dict.get_size() * 2 * sizeof(uint64_t)	/ (1ull << 20) << " MB)   "
		<< "patterns: " << v_kmer_patterns.size() - free_pattern_ids.size()
		<< " (" << tot_pat_size / (1ull << 20) << " MB)"
		<< "   -> "<< num_calc << " (" << num_calc / 1000000 << ")"
		<< " cv memory: " << cv_memory / 1000000 << "MB"
		<< " ht memory: " << ht_memory / 1000000 << "MB"
		<< " time: " << file_time.count()
		<< endl;
#else
	cout << "*** dict: " << hm_kmer_dict.get_size()
		<< " (" << hm_kmer_dict.get_size() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "patterns: " << v_kmer_patterns.size()
		<< " mem_pat: " << mem_pattern_desc / 1000000 << "MB"
		<< " mem_os_pat: " << mem_os_pattern_desc / 1000000 << "MB"
		<< " ht memory: " << ht_memory / 1000000 << "MB"
		<< " time: " << file_time.count()
		<< endl;
#endif

	fflush(stdout);
}

void store_pat_sizes(int num)
{
	if ((num < 500 && num % 10 != 0) || (num >= 500 && num % 100 != 0))
		return;

	FILE *stat;

	if (num == 0)
		stat = fopen("kmer_db.stat", "wb");
	else
		stat = fopen("kmer_db.stat", "ab");


	fprintf(stat, "\n************ %d \n", num);
	for (int i = 0; i < pat_sizes.size(); ++i)
	{
		auto x = pat_sizes[i];
		if(x)
			fprintf(stat, "    %5d:%7d\n", i, x);
	}

	fclose(stat);
}


int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "Usage: kmer_db <file_list>\n";
		cout << "   file_list - list with file_names to process\n";
		return 0;
	}

	load_file_list(string(argv[1]));

//	kmer_patterns[0] = make_pair(0, vid_t());
	hm_kmer_dict.set_special_keys((unsigned long long) - 1, (unsigned long long) - 2);
//	hm_current_mappings.set_special_keys((unsigned long long) - 1, (unsigned long long) - 2);
//	hm_kmer_patterns.set_special_keys((unsigned long long) - 1, (unsigned long long) - 2);
//	v_current_mappings.resize(16, 0);

//	hm_kmer_patterns.insert(0, make_pair(0, vid_t()));

//	v_kmer_patterns.resize(16, make_pair(0, vid_t()));
#ifdef PLAIN_VECTORS
	v_kmer_patterns.resize(16, pattern_t{ 0, vid_t(), 0 });
#else
	v_kmer_patterns.push_back(pattern_t{ 0, new vid_t()});
#endif
//	v_kmer_patterns[0] = make_pair(0, vid_t());
	for (uint32_t i = 1; i < v_kmer_patterns.size(); ++i)
		free_pattern_ids.push(i);

/*	v_patterns.resize(16);
	v_patterns[0] = pattern_t{ 0, 0, 0 };*/

	pat_sizes.resize(kmc_file_list.size() + 1, 0);

	for (size_t i = 0; i < kmc_file_list.size(); ++i)
	{
		cout << "Processing " << i << " of " << kmc_file_list.size() << "  -  " << kmc_file_list[i] << endl;

		auto start = std::chrono::high_resolution_clock::now();
		if (!process_single_file(i, kmc_file_list[i]))
		{
			cout << "Error\n";
			break;
		}
		auto end = std::chrono::high_resolution_clock::now();

		file_time = end - start;

		show_progress();
		store_pat_sizes(i);
	}

	store_pat_sizes(1000000);

#ifdef WIN32
//	getchar();
#endif

    return 0;
}

