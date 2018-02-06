
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

	bool useMinhash = false;
	it = find(params.begin(), params.end(), "-mh-input"); // minhash mode
	if (it != params.end()) {
		useMinhash = true;
		params.erase(it);
	}

	// all modes need at least 2 parameters
	if (params.size() >= 2) {
		const string& mode = params[0];
		
		if (params.size() == 3 && (mode == "build")) {
			cout << "Database building mode" << endl;
			return runBuildDatabase(params[1], params[2], useMinhash);
		}
		else if (params.size() == 3 && mode == "all2all") {
			cout << "All versus all comparison" << endl;
			return runAllVsAll(params[1], params[2]);
		}
		else if (params.size() == 4 && mode == "one2all") {
			cout << "One versus all comparison" << endl;
			return runOneVsAll(params[1], params[2], params[3], useMinhash);
		}
		else if (params.size() == 3 && mode == "list-patterns") {
			cout << "Listing all patterns" << endl;
			return runListPatterns(params[1], params[2]);
		}
		else if (params.size() == 3 && mode == "minhash" && std::atof(params[2].c_str()) > 0) {
			cout << "Min-hashing k-mers" << endl;
			return runMinHash(params[1], atof(params[2].c_str()));
		}
		else if (params.size() == 2 && mode == "distance") {
			cout << "Calculating distance measures" << endl;
			return runDistanceCalculation(params[1]);
		}
	}

	showInstructions();
	return 0;
}



int Console::runMinHash(const std::string& multipleKmcSamples, float fraction) {
	cout << "Minhashing samples..." << endl;

	std::chrono::duration<double> loadingTime, processingTime;

	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<MinHashFilter>(fraction, 0);

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
			file.store(task->filePath, *task->kmers, task->kmerLength);
		}

		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	return 0;
}



/****************************************************************************************************************************************************/

int Console::runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename, bool loadMinhash) {

	cout << "Processing samples..." << endl;
	
	LOG_DEBUG << "Creating FastKmerDb object" << endl;
	FastKmerDb* db = new FastKmerDb(numThreads);

	std::chrono::duration<double> loadingTime, processingTime;
	
	LOG_DEBUG << "Creating Loader object..." << endl;

	Loader loader(nullptr, loadMinhash, numThreads);
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
			db->addKmers(task->sampleName, *task->kmers, task->kmerLength);
			show_progress(*db);
		}
		
		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	loader.waitForPrefetch();

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
		<< "Number of k-mers: " << db->getKmersCount() << endl
		<< "K-mer length: " << db->getKmerLength() << endl;

	std::chrono::duration<double> dt;

	cout << "Serializing database...";
	auto start = std::chrono::high_resolution_clock::now();
	std::ofstream ofs;
	ofs.open(dbFilename, std::ios::binary);
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
		<< "Number of k-mers: " << db->getKmersCount() << endl
		<< "K-mer length: " << db->getKmerLength() << endl;


	cout << "Calculating matrix of common k-mers...";
	start = std::chrono::high_resolution_clock::now();
	LowerTriangularMatrix<uint32_t> matrix;
	db->calculateSimilarity(matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing matrix of common k-mers in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length, " << db->getKmerLength() << endl;
	std::copy(db->getSampleNames().cbegin(), db->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl <<  "total-kmers:" << endl;
	std::copy(db->getSampleKmersCount().cbegin(), db->getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl << "common-kmers:" << endl;
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

int Console::runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile, bool useMinhash) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads);

	std::chrono::duration<double> dt;
	cout << "Loading sample kmers...";

	auto start = std::chrono::high_resolution_clock::now();
	FastKmerDb sampleDb(numThreads);

	std::vector<kmer_t> kmers;
	uint32_t kmerLength;
	KmcFileWrapper file(nullptr);

	if (!file.open(singleKmcSample, useMinhash) || !file.load(kmers, kmerLength)) {
		cout << "FAILED";
		return -1;
	}
	sampleDb.addKmers(singleKmcSample, kmers, kmerLength);
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
		<< "Number of k-mers: " << db.getKmersCount() << endl
		<< "K-mer length: " << db.getKmerLength() << endl;


	cout << "Calculating similarity vector...";
	start = std::chrono::high_resolution_clock::now();
	std::vector<uint32_t> sims;
	db.calculateSimilarity(sampleDb, sims);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing similarity vector in " << similarityFile << "...";
	std::ofstream ofs(similarityFile);
	
	ofs << "sample-kmers, " << sampleDb.getSampleKmersCount().front() << endl
		<< "kmer-length, " << db.getKmerLength() << endl;
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl << "total-kmers:" << endl;
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl << "common-kmers:" << endl;
	std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));

	ofs.close();
	cout << "OK" << endl;

	return 0;
}


/****************************************************************************************************************************************************/

int Console::runDistanceCalculation(const std::string& similarityFilename) {

	std::vector<size_t> kmersCount;
	uint32_t kmerLength;
	uint64_t sampleKmersCount = 0;

	cout << "Loading file with common k-mer counts" << similarityFilename << "...";
	ifstream similarityFile(similarityFilename);
	if (!similarityFile) {
		cout << "FAILED" << endl;
		return -1;
	}
	cout << "OK" << endl;

	cout << "Calculating distances...";

	
	std::vector<std::string> metricNames = { "jaccard", "min", "max", "cosine", "mash" };
	std::vector<std::function<double(size_t, size_t, size_t, int)>> metrics = {
		[](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / (cnt1 + cnt2 - common); }, // jaccard
		[](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::min(cnt1,cnt2); }, // min
		[](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::max(cnt1,cnt2); }, // max
		[](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / sqrt(cnt1 * cnt2); }, // cosine
		[](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double {
			double d_jaccard = (double)common / (cnt1 + cnt2 - common); 
			return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1)); } // mash distance
	};

	std::vector<std::ofstream> files(metricNames.size());
	for (int i = 0; i < files.size(); ++i) {
		files[i].open(similarityFilename + "." + metricNames[i]);
	}


	string in;
	getline(similarityFile, in);  
	istringstream iss(in);

	// check if one vs all mode (number of sample kmers specified)
	iss >> in;
	if (in == "sample-kmers,") {
		iss >> sampleKmersCount;
		getline(similarityFile, in);
		iss = istringstream(in);
		iss >> in;
	}

	iss >> kmerLength; // get kmer length
	getline(similarityFile, in); // copy sample names to outout file
	for (auto & f : files) { f << in << endl; }
	getline(similarityFile, in); // ignore

	getline(similarityFile, in); // get number of kmers for all samples
	std::replace(in.begin(), in.end(), ',', ' ');
	
    istringstream iss2(in);
	std::copy(std::istream_iterator<size_t>(iss2), std::istream_iterator<size_t>(), std::back_inserter(kmersCount));
	getline(similarityFile, in); // ignore description row

	for (int i = 0; getline(similarityFile, in); ++i) {
		std::replace(in.begin(), in.end(), ',', ' ');
		istringstream iss(in);
		size_t intersection;

		uint64_t i_kmersCount = (sampleKmersCount > 0) ? sampleKmersCount : kmersCount[i];
			
		for (int j = 0; iss >> intersection; ++j) {
			// calculate all metrices
			for (int i = 0; i < metrics.size(); ++i) {
				files[i] << metrics[i](intersection, i_kmersCount, kmersCount[j], kmerLength) << ","; 
			}
		}
		
		for (auto & f : files) { f << endl; }
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
		<< "Number of k-mers: " << db.getKmersCount() << endl
		<< "K-mer length: " << db.getKmerLength() << endl;

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
		<< "\t   sample_list (input) - file containing list of k-mer files (raw or minhashed)" << endl
		<< "\t   database (output) - file with generated k-mer database" << endl

		<< "Minhashing k-mers:" << endl
		<< "\t kmer_db minhash <sample_list> <fraction>" << endl
		<< "\t   sample_list (input) - file containing list of k-mer files (raw)" << endl
		<< "\t   fraction (input) - fraction of kmers passing the filter" << endl

		<< "Calculating number of common k-mers for all the samples in the database:" << endl
		<< "\t kmer_db all2all <database> <common_matrix>" << endl
		<< "\t   database (input) - k-mer database file" << endl
		<< "\t   commony_matrix (output) - file containing matrix with numbers of common k-mers" << endl

		<< "Calculating number of common kmers between single sample and all the samples in the database:" << endl
		<< "\t kmer_db one2all <database> <sample> <common_vector>" << endl
		<< "\t   database (input) - k-mer database file" << endl
		<< "\t   sample (input) - k-mer file for a sample" << endl
		<< "\t   common_vector (output) - file containing vector with numbers of common k-mers" << endl

		<< "Calculating similarity/distance metrices (jaccard index and Mash distance) on the basis of common k-mers:" << endl
		<< "\t kmer_db distance <common_file>" << endl
		<< "\t   common_file (input) - file containing matrix/vector with numbers of common k-mers" << endl
		<< "This mode generates two files named <common_file>.jaccard and <common_file>.mash" << endl;
}

