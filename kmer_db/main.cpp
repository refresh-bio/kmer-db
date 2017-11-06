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
#include <iomanip>
#include <fstream>
#include <iterator>

using namespace std;

//#define COMPARE

#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif

#include "kmer_db.h"
#include "tests.h"

//map<uint32_t, uint32_t> pat_sizes;			// (rozmiar wzorca, liczba wyst wzorcow z takim rozmiarem)
vector<uint32_t> pat_sizes;						// rozmiar wzorca -> liczba wyst wzorcow z takim rozmiarem

												// Functions
bool load_file_list(string file_list);
void show_progress(const AbstractKmerDb &db);
void store_pat_sizes(int num);


std::string formatLargeNumber(uint64_t num) {
	std::string out = "";

	do {
		uint64_t part = num % 1000LL;
		num = num / 1000LL;

		if (num > 0) {
			ostringstream oss;
			oss << "," << setw(3) << setfill('0') << part;
			out = oss.str() + out;
		}
		else {
			out = std::to_string(part) + out;
		}

	} while (num > 0);

	return out;
}

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

	cout << "dict= " << formatLargeNumber(db.getKmersCount())
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "\t patterns= " << formatLargeNumber(db.getPatternsCount())
		<< "\t patterns mem= " << formatLargeNumber(db.getPatternBytes()) 
		<< "\t ht mem= " << formatLargeNumber(db.getHashtableBytes())
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
	
	pat_sizes.resize(kmc_file_list.size() + 1, 0);
	//kmc_file_list.resize(50);

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
					if (!fast_db.loadKmers(kmc_file_list[file_id], kmersCollections[tid])) {
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
				naive_db.addKmers(file_id, kmersCollections[tid]);
				dt = std::chrono::high_resolution_clock::now() - start;
				cout << "Naive: time=" << dt.count() << ", ";
				naiveTime += dt;
				show_progress(naive_db);
#endif

				start = std::chrono::high_resolution_clock::now();
				fast_db.addKmers(file_id, kmersCollections[tid]);
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


	cout << "DETAILED RESULTS:" << endl
		<< "Hashtable find: " << fast_db.hashtableFindTime.count() << endl
		<< "Hashatable add: " << fast_db.hashtableAddTime.count() << endl
		<< "Sort time: " << fast_db.sortTime.count() << endl
		<< "Pattern extension time: " << fast_db.extensionTime.count() << endl;


	std::ofstream ofs("kmer.db", std::ios::binary);
	fast_db.serialize(ofs);
	ofs.close();

	//Tests::testSerialization(fast_db);
	//Tests::testDistanceMatrix(fast_db, "d:/distances-fast.txt");


#ifdef COMPARE
	// compare naive and fast implementation
	cout << "NAIVE COMPARISON" << endl;

	Tests::comparePatterns(naive_db, fast_db, naive_db.getKmers());
	Tests::testDistanceMatrix(fast_db, "d:/distances-naive.txt");
	

#endif


	store_pat_sizes(1000000);

#ifdef WIN32
	//	getchar();
#endif

	return 0;
}

