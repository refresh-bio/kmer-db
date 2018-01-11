
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

	LOG_VERBOSE << "dict= " << Log::formatLargeNumber(db.getKmersCount())
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "\t patterns= " << Log::formatLargeNumber(db.getPatternsCount())
		<< "\t patterns mem= " << Log::formatLargeNumber(db.getPatternBytes())
		<< "\t ht mem= " << Log::formatLargeNumber(db.getHashtableBytes())
		<< endl;
}

/****************************************************************************************************************************************************/

int Console::parse(int argc, char** argv) {

	cout << "kmer-db version 1.0" << endl << endl;
	
	numThreads = 0;

	std::vector<string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* c)->string { return c; });

	// search for switches
	auto it = find(params.begin(), params.end(), "-v"); // verbose mode
	if (it != params.end()) {
		Log::getInstance(Log::LEVEL_VERBOSE).enable();
		params.erase(it);
	}

	it = find(params.begin(), std::prev(params.end()), "-t"); // number of threads
	if (it != std::prev(params.end())) {
		numThreads = atof(std::next(it)->c_str());
		params.erase(it, it + 2);
	}

	// all modes need at least 3 parameters
	if (params.size() >= 3) {
		const string& mode = params[0];
		
		if (params.size() == 3 && mode == "build") {
			cout << "Database building mode" << endl;
			return runBuildDatabase(params[1], params[2]);
		}
		else if (params.size() == 3 && mode == "all2all") {
			cout << "All versus all comparison" << endl;
			return runAllVsAll(params[1], params[2]);
		}
		else if (params.size() == 4 && mode == "one2all") {
			cout << "One versus all comparison" << endl;
			return runOneVsAll(params[1], params[2], params[3]);
		}
		else if (params.size() == 3 && mode == "list-patterns") {
			cout << "Listing all patterns" << endl;
			return runListPatterns(params[1], params[2]);
		}
		else if (params.size() == 3 && mode == "minhash" && std::atof(params[2].c_str()) > 0) {
			cout << "Min-hashing k-mers" << endl;
			return runMinHash(params[1], atof(params[2].c_str()));
		}
		else if (params.size() == 3 && mode == "distance") {
			cout << "Calculating distance measures" << endl;
			return runDistanceCalculation(params[1], params[2]);
		}
	}

	showInstructions();
	return 0;
}



int Console::runMinHash(const std::string& multipleKmcSamples, float fraction) {
	cout << "Minhashing samples..." << endl;

	std::chrono::duration<double> loadingTime, processingTime;

	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<MinHashFilter>(fraction, 20);

	Loader loader(filter, false, numThreads);
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
			KmcFileWrapper file(nullptr); 
			file.store(task->filePath, *task->kmers);
		}

		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	return 0;
}



/****************************************************************************************************************************************************/

int Console::runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename) {

	cout << "Processing samples..." << endl;
	
	LOG_DEBUG << "Creating FastKmerDb object" << endl;
	FastKmerDb* db = new FastKmerDb(numThreads);

	std::chrono::duration<double> loadingTime, processingTime;
	
	LOG_DEBUG << "Creating Loader object..." << endl;

	Loader loader(nullptr, true, numThreads);
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
		}
		
		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	cout << endl << endl << "EXECUTION TIMES" << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Processing time: " << processingTime.count() << endl
		<< "\tHashatable resizing (serial): " << db->hashtableResizeTime.count() << endl
		<< "\tHashtable searching (parallel): " << db->hashtableFindTime.count() << endl
		<< "\tHashatable insertion (serial): " << db->hashtableAddTime.count() << endl
		<< "\tSort time (parallel): " << db->sortTime.count() << endl
		<< "\tPattern extension time (serial): " << db->extensionTime.count() << endl << endl
		<< "STATISTICS" << endl
		<< "Number of samples: " << db->getSamplesCount() << endl
		<< "Number of patterns: " << db->getPatternsCount() << endl
		<< "Number of k-mers: " << db->getKmersCount() << endl;

	std::chrono::duration<double> dt;

	cout << "Serializing database...";
	auto start = std::chrono::high_resolution_clock::now();
	std::ofstream ofs;
	ofs.open(dbFilename, std::ios::binary);
	db->serialize(ofs);
	ofs.close();

	ofs.open(dbFilename + ".log");
	for (int i = 0; i < db->getSamplesCount(); ++i) {
		ofs << db->getSampleNames()[i] << " " << db->getSampleKmersCount()[i] << endl;
	}
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
	FastKmerDb* db = new FastKmerDb(numThreads);;

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

	return 0;
}


/****************************************************************************************************************************************************/

int Console::runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads);

	std::chrono::duration<double> dt;
	cout << "Loading sample kmers...";

	auto start = std::chrono::high_resolution_clock::now();
	FastKmerDb sampleDb(numThreads);

	std::vector<kmer_t> kmers;
	KmcFileWrapper file(nullptr);

	if (!file.open(singleKmcSample, true) || !file.load(kmers)) {
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

	return 0;
}


/****************************************************************************************************************************************************/

int Console::runDistanceCalculation(const std::string& dbFilename, const std::string& similarityFilename) {
	
	std::string rerportFilename = dbFilename + ".log";
	
	std::ifstream report(rerportFilename);

	std::vector<size_t> kmersCount;
	
	cout << "Loading k-mer statistics from report " << rerportFilename << "...";
	if (!report) {
		cout << "FAILED" << endl
		 << "Loading k-mer statistics from database " << dbFilename << "...";
		std::ifstream dbFile(dbFilename, std::ios::binary);
		FastKmerDb db;
		 if (!dbFile || !db.deserialize(dbFile)) {
			 cout << "FAILED" << endl;
			 return -1;
		 }
		 else {
			 cout << "OK" << endl;
			 kmersCount = db.getSampleKmersCount();
		 }
	} else {
		cout << "OK" << endl;
		string name;
		size_t count;
		while (report >> name >> count) {
			kmersCount.push_back(count);
		}
	}


	cout << "Loading raw similarity file " << similarityFilename << "...";
	ifstream similarityFile(similarityFilename);
	if (!similarityFile) {
		cout << "FAILED" << endl;
		return -1;
	}
	cout << "OK" << endl;

	cout << "Calculating distances...";
	ofstream jaccardFile(similarityFilename + ".jaccard");
	ofstream mashFile(similarityFilename + ".mash");

	string in;
	getline(similarityFile, in); // ignore first row

	for (int i = 0; i < kmersCount.size() && getline(similarityFile, in); ++i) {
		std::replace(in.begin(), in.end(), ',', ' ');
		istringstream iss(in);
		size_t intersection;

		int j = 0;
		for (; j < i && iss >> intersection; ++j) {
			double d_jaccard = (double)intersection / (kmersCount[i] + kmersCount[j] - intersection);
			// fixme: fixed kmer length
			double d_mash = (d_jaccard == 0) ? 1.0 : (-1.0 / 18) * log((2 * d_jaccard) / (d_jaccard + 1)); 
			
			jaccardFile << d_jaccard << ",";
			mashFile << d_mash << ",";
		}
		
		for (; j < kmersCount.size(); ++j) {
			jaccardFile << ",0";
			mashFile << ",0";
		}
		jaccardFile << endl;
		mashFile << endl;
		++i;
	}

	cout << "OK" << endl;
}


/****************************************************************************************************************************************************/

int Console::runListPatterns(const std::string& dbFilename, const std::string& patternFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads);

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl
		<< "Number of samples: " << db.getSamplesCount() << endl
		<< "Number of patterns: " << db.getPatternsCount() << endl
		<< "Number of k-mers: " << db.getKmersCount() << endl;

	cout << "Storing patterns in " << patternFile << "...";
	std::ofstream ofs(patternFile);
	db.savePatterns(ofs);
	ofs.close();
	cout << "OK" << endl; 

	return 0;
}


/****************************************************************************************************************************************************/
void Console::showInstructions() {
	cout	<< "USAGE" << endl

		<< "Building k-mer database:" << endl
		<< "\t kmer_db build <sample_list> <database>" << endl
		<< "\t   sample_list (input) - file containing list of k-mer files (raw or min-hashed)" << endl
		<< "\t   database (output) - file with generated k-mer database" << endl

		<< "Min-hashing k-mers:" << endl
		<< "\t kmer_db minhash <sample_list> <fraction>" << endl
		<< "\t   sample_list (input) - file containing list of k-mer files (raw)" << endl
		<< "\t   fraction (input) - fraction of kmers passing the filter" << endl

		<< "Calculating similarity matrix for all the samples in the database:" << endl
		<< "\t kmer_db all2all <database> <similarity_matrix>" << endl
		<< "\t   database (input) - k-mer database file" << endl
		<< "\t   similarity_matrix (output) - file with similarity matrix" << endl

		<< "Calculating similarity of a new sample to all the samples in the database:" << endl
		<< "\t kmer_db one2all <database> <sample> <similarity_vector>" << endl
		<< "\t   database (input) - k-mer database file" << endl
		<< "\t   sample (input) - k-mer file for a sample" << endl
		<< "\t   similarity_matrix (output) - file with similarity matrix" << endl;
}

