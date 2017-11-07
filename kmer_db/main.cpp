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
#include "log.h"

		

void showInstructions() {
	cout << "kmer-db version 1.0" << endl
		<< "USAGE" << endl

		<< "Building k-mer database:" << endl
		<< "\t kmer_db --build <sample_list> <database>" << endl
		<< "\t   sample_list - directory with k-mer files or file containing list of k-mer files (input)" << endl
		<< "\t   database - file with generated k-mer database (output)" << endl

		<< "Calculating similarity matrix for all the samples in the database:" << endl
		<< "\t kmer_db --all2all <database> <similarity_matrix>" << endl
		<< "\t   database - k-mer database file (input)" << endl
		<< "\t   similarity_matrix - file with similarity matrix (output)" << endl

		<< "Calculating similarity of a new sample to all the samples in the database:" << endl
		<< "\t kmer_db --one2all <sample> <database> <similarity_vector>" << endl
		<< "\t   sample - k-mer file for a sample (input)" << endl
		<< "\t   database - k-mer database file (input)" << endl
		<< "\t   similarity_matrix - file with similarity matrix (output)" << endl;
}

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



int main(int argc, char **argv)
{
	FastKmerDb db;
	Log::getInstance(Log::LEVEL_DEBUG).enable();

	if (argc == 4 && string(argv[1]) == "--build") {
		
		vector<string> kmc_file_list;
		load_file_list(string(argv[2]), kmc_file_list);
		std::ofstream ofs(argv[3], std::ios::binary);

		kmc_file_list.resize(100);

		std::chrono::duration<double> loadingTime, fastTime;
		cout << "PROCESSING SAMPLES..." << endl;

		int n_threads = std::thread::hardware_concurrency();
		std::vector<std::vector<uint64_t>> kmersCollections(n_threads);
		std::vector<std::thread> readers(n_threads);

		for (size_t i = 0; i < kmc_file_list.size(); i += n_threads)
		{
			std::chrono::duration<double> dt;
			auto start = std::chrono::high_resolution_clock::now();

			for (int tid = 0; tid < n_threads; ++tid) {
				readers[tid] = std::thread([tid, i, n_threads, &kmersCollections, &kmc_file_list, &db]() {
					size_t file_id = i + tid;

					if (file_id < kmc_file_list.size()) {
						std::ostringstream oss;
						oss << kmc_file_list[file_id] << " (" << file_id + 1 << "/" << kmc_file_list.size() << ")...";
						if (!db.loadKmers(kmc_file_list[file_id], kmersCollections[tid])) {
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
					start = std::chrono::high_resolution_clock::now();
					db.addKmers(file_id, kmersCollections[tid]);
					dt = std::chrono::high_resolution_clock::now() - start;
					cout << "Fast: time=" << dt.count() << ", ";
					fastTime += dt;
					show_progress(db);
					cout << endl;
				}
			}
		}

		cout << endl << "EXECUTION TIMES" << endl
			<< "Loading k-mers: " << loadingTime.count() << endl
			<< "Total processing time: " << fastTime.count() << endl
			<< "\tHashatable resizing (serial): " << db.hashtableResizeTime.count() << endl
			<< "\tHashtable searching (parallel): " << db.hashtableFindTime.count() << endl
			<< "\tHashatable insertion (serial): " << db.hashtableAddTime.count() << endl
			<< "\tSort time (parallel): " << db.sortTime.count() << endl
			<< "\tPattern extension time (serial): " << db.extensionTime.count() << endl;

		db.serialize(ofs);
		ofs.close();
	}
	else if (argc == 4 && string(argv[1]) == "--all2all") {
		std::ifstream dbFile(argv[2], std::ios::binary);
		std::ofstream matrixFile(argv[3]);

		if (dbFile && matrixFile) {
			Array<uint32_t> matrix;
			db.deserialize(dbFile);
			db.calculateSimilarityMatrix(matrix);

		}

	}
	else if (argc == 5 && string(argv[1]) == "--one2all") {
	
	}
	else {
		showInstructions();
	}

	

	return 0;
}

