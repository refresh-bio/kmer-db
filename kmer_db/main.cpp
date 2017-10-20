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

#include "kmer_db.h"

using namespace std;





//map<uint32_t, uint32_t> pat_sizes;			// (rozmiar wzorca, liczba wyst wzorcow z takim rozmiarem)
vector<uint32_t> pat_sizes;						// rozmiar wzorca -> liczba wyst wzorcow z takim rozmiarem

												// Functions
bool load_file_list(string file_list);
void show_progress(const FastKmerDb &db);
void store_pat_sizes(int num);

// Wczytuje plik zawierajacy liste baz KMC do przetworzenia
bool load_file_list(const string& file_list, vector<string>& kmc_file_list)
{
	char s[1024];

	kmc_file_list.clear();

	FILE *f = fopen(file_list.c_str(), "rt");
	if (!f)
		return false;

	while (!feof(f))
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


// Pokazuje stan
void show_progress(const FastKmerDb &db)
{
	size_t tot_pat_size = 0;
	size_t num_calc = 0;			// Liczba operacji przy wyznaczaniu macierzy podobieñstwa

	cout << "*** dict: " << db.getKmersCount()
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "patterns: " << db.getPatternsCount()
		<< " mem_pat: " << mem_pattern_desc / 1000000 << "MB"
		<< " mem_os_pat: " << mem_os_pattern_desc / 1000000 << "MB"
		<< " ht memory: " << ht_memory / 1000000 << "MB"
		<< endl;

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
		if (x)
			fprintf(stat, "    %5d:%7d\n", i, x);
	}

	fclose(stat);
}


int main(int argc, char **argv)
{
	vector<string> kmc_file_list;
	
	if (argc < 2)
	{
		cout << "Usage: kmer_db <file_list>\n";
		cout << "   file_list - list with file_names to process\n";
		return 0;
	}

	load_file_list(string(argv[1]), kmc_file_list);

	FastKmerDb fast_db;
	NaiveKmerDb naive_db;
	std::vector<uint64_t> kmers;

	pat_sizes.resize(kmc_file_list.size() + 1, 0);
	kmc_file_list.resize(2);

	std::chrono::duration<double> loadingTime, naiveTime, fastTime;

	cout << "PROCESSING SAMPLES..." << endl;

	for (size_t i = 0; i < kmc_file_list.size(); ++i)
	{
		cout << kmc_file_list[i] << "(" << i + 1 << "/" << kmc_file_list.size() << "): ";

		std::chrono::duration<double> dt;

		auto start = std::chrono::high_resolution_clock::now();
		if (!fast_db.extractKmers(kmc_file_list[i], kmers)) {
			cout << "Error processing sample" << endl;
			break;
		}
		dt = std::chrono::high_resolution_clock::now() - start;
		cout << "loading time=" << dt.count();
		loadingTime += dt;

		start = std::chrono::high_resolution_clock::now();
		naive_db.addKmers(i, kmers);
		dt = std::chrono::high_resolution_clock::now() - start;
		cout << ", naive time=" << dt.count();
		naiveTime += dt;

		start = std::chrono::high_resolution_clock::now();
		fast_db.addKmers(i, kmers);
		dt = std::chrono::high_resolution_clock::now() - start;
		cout << ", fast time=" << dt.count();
		fastTime += dt;
		
		cout << endl;
		show_progress(fast_db);
		store_pat_sizes(i);
	}

	cout << endl << "EXECUTION TIMES" << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Naive processing: " << naiveTime.count() << endl
		<< "Fast processing: " << fastTime.count() << endl;

	// compare naive and fast implementation
	cout << "Comparing naive and fast implementations" << endl;
	std::vector<sample_id_t> samples;
	for (auto it = naive_db.kmers2samples.begin(); it != naive_db.kmers2samples.end(); ++it) {
		auto& naiveSamples = it->second;
		fast_db.getKmerSamples(it->first, samples);

		auto mis = std::mismatch(naiveSamples.begin(), naiveSamples.end(), samples.begin());

		if (mis.first != naiveSamples.end() || mis.second != samples.end()) {
			cout << "k-mer: " << it->first << endl;
		}

	}
	cout << "done" << endl;

	store_pat_sizes(1000000);

#ifdef WIN32
	//	getchar();
#endif

	return 0;
}

