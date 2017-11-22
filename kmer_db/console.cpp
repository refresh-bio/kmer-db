
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
#include "loader.h"


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

	cout << "Processing samples..." << endl;
	
	LOG_DEBUG << "Creating FastKmerDb object" << endl;
	FastKmerDb* db = new FastKmerDb();
	std::ofstream ofs(dbFilename, std::ios::binary);

	std::chrono::duration<double> loadingTime, processingTime;
	
	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<NullFilter>();
	//auto filter = std::make_shared<MinHashFilter>(0.1, 20);

	Loader loader(filter);
	loader.configure(multipleKmcSamples);

	loader.initPrefetch();

	LOG_DEBUG << "Starting loop..." << endl;

	for (;;) {
		auto start = std::chrono::high_resolution_clock::now();
		loader.waitForPrefetch();
		loader.initLoad();
		loader.waitForLoad();
		loadingTime += std::chrono::high_resolution_clock::now() - start;
		
		start = std::chrono::high_resolution_clock::now();
		loader.initPrefetch();
		if (!loader.getLoadedTasks().size()) {
			break;
		}
		
		for (const auto& entry : loader.getLoadedTasks()) {
			auto task = entry.second;
			db->addKmers(task->sampleName, *task->kmers);
			show_progress(*db);
			cout << endl;
		}
		
		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	cout << endl << "EXECUTION TIMES" << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Processing time: " << processingTime.count() << endl
		<< "\tHashatable resizing (serial): " << db->hashtableResizeTime.count() << endl
		<< "\tHashtable searching (parallel): " << db->hashtableFindTime.count() << endl
		<< "\tHashatable insertion (serial): " << db->hashtableAddTime.count() << endl
		<< "\tSort time (parallel): " << db->sortTime.count() << endl
		<< "\tPattern extension time (serial): " << db->extensionTime.count() << endl;

	std::chrono::duration<double> dt;

	cout << "Serializing database...";
	auto start = std::chrono::high_resolution_clock::now();
	db->serialize(ofs);
	ofs.close();
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Releasing memory...";
	start = std::chrono::high_resolution_clock::now();
	delete db;
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	return 0;
}

/****************************************************************************************************************************************************/

int Console::runAllVsAll(const std::string& dbFilename, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	std::ofstream ofs(similarityFile);
	FastKmerDb* db = new FastKmerDb();;

	std::chrono::duration<double> dt;
	cout << "Loading k-mer database " << dbFilename << "...";
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db->deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl
		<< "Number of samples: " << db->getSamplesCount() << endl
		<< "Number of patterns: " << db->getPatternsCount() << endl
		<< "Number of k-mers: " << db->getKmersCount() << endl;


	cout << "Calculating similarity matrix...";
	start = std::chrono::high_resolution_clock::now();
	LowerTriangularMatrix<uint32_t> matrix;
	db->calculateSimilarity(matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing similarity matrix in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	std::copy(db->getSampleNames().begin(), db->getSampleNames().end(), ostream_iterator<string>(ofs, ","));
	ofs << endl;
	matrix.save(ofs);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Releasing memory...";
	start = std::chrono::high_resolution_clock::now();
	delete db;
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;
}


/****************************************************************************************************************************************************/

int Console::runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db;

	std::chrono::duration<double> dt;
	cout << "Loading sample kmers...";

	auto start = std::chrono::high_resolution_clock::now();
	FastKmerDb sampleDb;
	Loader loader(make_shared<NullFilter>());
	std::vector<kmer_t> kmers;
	if (!loader.loadKmers(singleKmcSample, kmers)) {
		cout << "FAILED";
		return -1;
	}
	sampleDb.addKmers(singleKmcSample, kmers);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl
		<< "Number of k-mers: " << sampleDb.getKmersCount() << endl;

	cout << "Loading k-mer database " << dbFilename << "...";
	start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl
		<< "Number of samples: " << db.getSamplesCount() << endl
		<< "Number of patterns: " << db.getPatternsCount() << endl
		<< "Number of k-mers: " << db.getKmersCount() << endl;


	cout << "Calculating similarity vector...";
	start = std::chrono::high_resolution_clock::now();
	std::vector<uint32_t> sims;
	db.calculateSimilarity(sampleDb, sims);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

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

	return 0;
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
		<< "\t kmer_db --one2all <database> <sample> <similarity_vector>" << endl
		<< "\t   database - k-mer database file (input)" << endl
		<< "\t   sample - k-mer file for a sample (input)" << endl
		<< "\t   similarity_matrix - file with similarity matrix (output)" << endl;
}

