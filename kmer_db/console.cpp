
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iterator>

#include "console.h"
#include "kmer_db.h"

using namespace std;


// Pokazuje stan
void show_progress(const AbstractKmerDb &db)
{
	size_t tot_pat_size = 0;
	size_t num_calc = 0;			// Liczba operacji przy wyznaczaniu macierzy podobieñstwa

	cout << "dict= " << Log::formatLargeNumber(db.getKmersCount())
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "\t patterns= " << Log::formatLargeNumber(db.getPatternsCount())
		<< "\t patterns mem= " << Log::formatLargeNumber(db.getPatternBytes())
		<< "\t ht mem= " << Log::formatLargeNumber(db.getHashtableBytes())
		<< endl;

	fflush(stdout);
}

/****************************************************************************************************************************************************/

int Console::parse(int argc, char** argv) {

	cout << "kmer-db version 1.0" << endl;
	
		if (argc == 4 && string(argv[1]) == "--build") {
		cout << "Database building mode" << endl;
		return runBuildDatabase(argv[2], argv[3]);
	}
	else if (argc == 4 && string(argv[1]) == "--all2all") {
		cout << "All versus all comparison" << endl;
		return runAllVsAll(argv[2], argv[3]);
	}
	else if (argc == 5 && string(argv[1]) == "--one2all") {
		cout << "One versus all comparison" << endl;
		return runOneVsAll(argv[2], argv[3], argv[4]);
	}
	else if (argc == 4 && string(argv[1]) == "--list-patterns") {
		cout << "Listing all patterns" << endl;
		return runListPatterns(argv[2], argv[3]);
	}
	else {
		showInstructions();
		return 0;
	}
}


/****************************************************************************************************************************************************/

int Console::runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename) {

	FastKmerDb db;
	vector<string> kmc_file_list;
	loadFileList(multipleKmcSamples, kmc_file_list);
	std::ofstream ofs(dbFilename, std::ios::binary);

	//		kmc_file_list.resize(100);

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
				string sample = kmc_file_list[file_id];

				size_t pos = sample.find_last_of("/\\");
				if (pos != string::npos) {
					sample = sample.substr(pos + 1);
				}

				start = std::chrono::high_resolution_clock::now();
				db.addKmers(sample, kmersCollections[tid]);
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

	return 0;
}

/****************************************************************************************************************************************************/

int Console::runAllVsAll(const std::string& dbFilename, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	std::ofstream ofs(similarityFile);
	FastKmerDb db;

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl
		<< "Number of samples: " << db.getSamplesCount() << endl
		<< "Number of patterns: " << db.getPatternsCount() << endl
		<< "Number of k-mers: " << db.getKmersCount() << endl;

	cout << "Calculating similarity matrix...";
	Array<uint32_t> matrix;
	db.calculateSimilarity(matrix);
	cout << "OK" << endl;

	cout << "Storing similarity matrix in " << similarityFile << "...";
	std::copy(db.getSampleNames().begin(), db.getSampleNames().end(), ostream_iterator<string>(ofs, ","));
	ofs << endl;
	
	matrix.save(ofs);
	cout << "OK" << endl;
}


/****************************************************************************************************************************************************/

int Console::runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db;


	cout << "Loading sample kmers...";
	FastKmerDb sampleDb;
	std::vector<kmer_t> kmers;
	if (!sampleDb.loadKmers(singleKmcSample, kmers)) {
		cout << "FAILED";
		return -1;
	}
	sampleDb.addKmers(singleKmcSample, kmers);
	cout << "OK" << endl
		<< "Number of k-mers: " << sampleDb.getKmersCount() << endl;

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl
		<< "Number of samples: " << db.getSamplesCount() << endl
		<< "Number of patterns: " << db.getPatternsCount() << endl
		<< "Number of k-mers: " << db.getKmersCount() << endl;


	cout << "Calculating similarity vector...";
	std::vector<uint32_t> sims;
	db.calculateSimilarity(sampleDb, sims);
	cout << "OK" << endl;

	cout << "Storing similarity vector in " << similarityFile << "...";
	std::ofstream ofs(similarityFile);
	
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl;
	std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));

	ofs.close();
	cout << "OK" << endl;
}

/****************************************************************************************************************************************************/

int Console::runListPatterns(const std::string& dbFilename, const std::string& patternFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db;

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl
		<< "Number of samples: " << db.getSamplesCount() << endl
		<< "Number of patterns: " << db.getPatternsCount() << endl
		<< "Number of k-mers: " << db.getKmersCount() << endl;

	cout << "Storing patterns in in " << patternFile << "...";
	std::ofstream ofs(patternFile);
	db.savePatterns(ofs);
	ofs.close();
}

/****************************************************************************************************************************************************/

bool Console::loadFileList(const std::string& multipleKmcSamples, std::vector<std::string>& kmcFileList) {

	std::ifstream ifs(multipleKmcSamples);

	string fname;
	while (ifs >> fname) {
		kmcFileList.push_back(fname);
	}

	sort(kmcFileList.begin(), kmcFileList.end());
	kmcFileList.erase(unique(kmcFileList.begin(), kmcFileList.end()), kmcFileList.end());

	return !kmcFileList.empty();
}

/****************************************************************************************************************************************************/
void Console::showInstructions() {
	cout	<< "USAGE" << endl

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

