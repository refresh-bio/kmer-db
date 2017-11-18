#include "tests.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "console.h"

using namespace std;

bool load_file_list(const string& file_list, vector<string>& kmc_file_list);
void show_progress(const AbstractKmerDb &db);

/*
void Tests::compareWithNaive(const string& fname) {
	vector<string> kmc_file_list;

	Console console;
	console.loadFileList(fname, kmc_file_list);

	FastKmerDb fast_db;
	NaiveKmerDb naive_db;
	

	//kmc_file_list.resize(100);

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

				string sample = kmc_file_list[file_id];

				size_t pos = sample.find_last_of("/\\");
				if (pos != string::npos) {
					sample = sample.substr(pos + 1);
				}

				start = std::chrono::high_resolution_clock::now();
				naive_db.addKmers(sample, kmersCollections[tid]);
				dt = std::chrono::high_resolution_clock::now() - start;
				cout << "Naive: time=" << dt.count() << ", ";
				naiveTime += dt;
				show_progress(naive_db);

				start = std::chrono::high_resolution_clock::now();
				fast_db.addKmers(sample, kmersCollections[tid]);
				dt = std::chrono::high_resolution_clock::now() - start;
				cout << "Fast: time=" << dt.count() << ", ";
				fastTime += dt;
				show_progress(fast_db);

				cout << endl;
			}
		}
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


	// compare naive and fast implementation
	cout << "NAIVE COMPARISON" << endl;

	Tests::comparePatterns(naive_db, fast_db, naive_db.getKmers());
	Tests::testDistanceMatrix(fast_db, "d:/distances-naive.txt");

}

*/


void Tests::comparePatterns(const AbstractKmerDb& db1, const AbstractKmerDb& db2, const std::vector<kmer_t>& kmers) {

	cout << "Comparing databases..." << endl;

	int i = 0;

	std::vector<sample_id_t> samples1;
	std::vector<sample_id_t> samples2;

	for (auto kmer : kmers) {
		db1.mapKmers2Samples(kmer, samples1);
		db2.mapKmers2Samples(kmer, samples2);
			
		bool eq = (samples1 == samples2);

		if (!eq) {
			std::cout << "k-mer: " << kmer << std::endl;
		}

	}
}


void Tests::testSerialization(const FastKmerDb& db) {
	
	cout << "Serializing..." << endl;

	std::ofstream ofs("D:/kmer-db.bin", std::ios::binary);
	db.serialize(ofs);
	ofs.close();

	cout << "Deserializing..." << endl;

	std::ifstream ifs("D:/kmer-db.bin", std::ios::binary);
	FastKmerDb copy;
	copy.deserialize(ifs);
	ifs.close();

	cout << "Comparing...";
	comparePatterns(db, copy, db.getKmers());

	cout << "Done!" << endl;

}


 void Tests::testDistanceMatrix(const AbstractKmerDb& db, const string& fname) {

	 cout << "Calculating distance matrix..." << endl;

	 std::ofstream file(fname);

	 file << endl << "Distance matrix:" << endl;
	 Array<uint32_t> matrix;
//	 db.calculateSimilarity(matrix);
	 cout << "Saving distance matrix...";
	 for (int i = 0; i < matrix.size(); ++i) {
		 for (int j = 0; j < matrix.size(); ++j) {
			 file << setw(10) << matrix[i][j];
		 }
		 file << endl;
	 }

	 cout << "done!" << endl;

/*	 fileNaive << endl << "Histogram: " << endl;
	 auto stats = naive_db.getPatternsStatistics();
	 for (auto s : stats) {
		 fileNaive << setw(10) << s.second << ": ";
		 copy(s.first.begin(), s.first.end(), ostream_iterator<sample_id_t>(fileNaive, ","));
		 fileNaive << endl;
	 }
	 */
}