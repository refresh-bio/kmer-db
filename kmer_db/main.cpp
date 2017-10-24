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
#include <sstream>

using namespace std;

//#define COMPARE

#ifndef WIN32
#include <parallel/algorithm>
#endif

#include "kmer_db.h"

//map<uint32_t, uint32_t> pat_sizes;			// (rozmiar wzorca, liczba wyst wzorcow z takim rozmiarem)
vector<uint32_t> pat_sizes;						// rozmiar wzorca -> liczba wyst wzorcow z takim rozmiarem

												// Functions
bool load_file_list(string file_list);
void show_progress(const AbstractKmerDb &db);
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
void show_progress(const AbstractKmerDb &db)
{
	size_t tot_pat_size = 0;
	size_t num_calc = 0;			// Liczba operacji przy wyznaczaniu macierzy podobieñstwa

	cout << "dict= " << db.getKmersCount()
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "patterns= " << db.getPatternsCount()
		<< " patterns mem= " << db.getPatternMem() / 1000000 << "MB"
		<< " ht mem= " << db.getHashtableMem() / 1000000 << "MB"
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
	//std::vector<uint64_t> kmers;

	pat_sizes.resize(kmc_file_list.size() + 1, 0);
	kmc_file_list.resize(50);

	std::chrono::duration<double> loadingTime, naiveTime, fastTime;

	cout << "PROCESSING SAMPLES..." << endl;
	
	int n_threads = std::thread::hardware_concurrency();
	std::vector<std::vector<uint64_t>> kmersCollections(n_threads);
	std::vector<std::thread> readers(n_threads);


	for (size_t i = 0; i < kmc_file_list.size(); i += n_threads)
	{
		std::chrono::duration<double> dt;
		auto start = std::chrono::high_resolution_clock::now();
		
		for (int tid = 0; tid < n_threads; ++tid) {
			readers[tid] = std::thread([tid, i, n_threads, &kmersCollections, &kmc_file_list, &fast_db, &naive_db]() {
				size_t file_id = i + tid;

				if (file_id < kmc_file_list.size()) {
					std::ostringstream oss;
					oss << kmc_file_list[file_id] << " (" << file_id + 1 << "/" << kmc_file_list.size() << ")...";
					if (!fast_db.extractKmers(kmc_file_list[file_id], kmersCollections[tid])) {
						oss << "Error processing sample" << endl;
					}
					else {
						oss << "done!" << endl;
					}
					cout << oss.str();
				}
			});
		}
		

		for (int tid = 0; tid < n_threads; ++tid) {
			readers[tid].join();
		}
	
		dt = std::chrono::high_resolution_clock::now() - start;
		cout << "Processed " << n_threads << " files in: " << dt.count() << endl;
		loadingTime += dt;
		
		for (int tid = 0; tid < n_threads; ++tid) {
			size_t file_id = i + tid;

			if (file_id < kmc_file_list.size()) {
#ifdef COMPARE
				start = std::chrono::high_resolution_clock::now();
				naive_db.addKmers(i, kmersCollections[tid]);
				dt = std::chrono::high_resolution_clock::now() - start;
				cout << "Naive: time=" << dt.count() << ", ";
				naiveTime += dt;
				show_progress(naive_db);
#endif

				start = std::chrono::high_resolution_clock::now();
				fast_db.addKmers(i, kmersCollections[tid]);
				dt = std::chrono::high_resolution_clock::now() - start;
				cout << "Fast: time=" << dt.count() << ", ";
				fastTime += dt;
				show_progress(fast_db);
				
				cout << endl;
			}
		}

		store_pat_sizes(i);
	}

	cout << endl << "EXECUTION TIMES" << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Naive processing: " << naiveTime.count() << endl
		<< "Fast processing: " << fastTime.count() << endl;

	// compare naive and fast implementation

#ifdef COMPARE
	cout << "Comparing naive and fast implementations" << endl;
	std::vector<sample_id_t> samples;
	int i = 0;

	for (auto it = naive_db.kmers2patternIds.begin(); it < naive_db.kmers2patternIds.end(); ++it, ++i) {
		
		if (naive_db.kmers2patternIds.is_free(*it)) {
			continue;
		}
		
		auto patternId = it->val;
		auto& naiveSamples = naive_db.patterns[patternId];
		
		fast_db.mapKmers2Samples(it->key, samples);

		bool eq = std::equal(naiveSamples.begin(), naiveSamples.end(), samples.begin(), samples.end());

		if (!eq) {
			cout << "k-mer: " << it->key << endl;

		}

	}
	cout << "done" << endl;
#endif

	store_pat_sizes(1000000);

#ifdef WIN32
	//	getchar();
#endif

	return 0;
}

