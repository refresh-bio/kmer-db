/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

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
#include <functional>

#include "console.h"
#include "kmer_db.h"
#include "loader.h"
#include "version.h"
#include "analyzer.h"


using namespace std;


// *****************************************************************************************
//
void show_progress(const FastKmerDb &db)
{
	size_t tot_pat_size = 0;
	size_t num_calc = 0;			

	LOG_VERBOSE << "dict= " << Log::formatLargeNumber(db.getKmersCount())
		<< " (" << db.getKmersCount() * 2 * sizeof(uint64_t) / (1ull << 20) << " MB)   "
		<< "\t patterns= " << Log::formatLargeNumber(db.getPatternsCount())
		<< "\t patterns mem= " << Log::formatLargeNumber(db.getPatternBytes())
		<< "\t ht mem= " << Log::formatLargeNumber(db.getHashtableBytes())
		<< endl;
}

// *****************************************************************************************
//
int Console::parse(int argc, char** argv) {

	cout << "Kmer-db version " << VERSION << " (" << DATE << ")" << endl 
		<< "S. Deorowicz, A. Gudys, M. Dlugosz, M. Kokot, and A. Danek (c) 2018" << endl << endl;
	
	numThreads = 0;
	cacheBufferMb = 8;

	std::vector<string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* c)->string { return c; });

	if (params.size()) {

		// search for switches
		auto it = find(params.begin(), params.end(), "-v"); // verbose mode
		if (it != params.end()) {
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
			params.erase(it);
		}

		it = find(params.begin(), std::prev(params.end()), "-t"); // number of threads
		if (it != std::prev(params.end())) {
			numThreads = (int) (atof(std::next(it)->c_str()));
			params.erase(it, it + 2);
		}

		it = find(params.begin(), std::prev(params.end()), "-buffer"); // size of temporary buffer in megabytes
		if (it != std::prev(params.end())) {
			cacheBufferMb = atoi(std::next(it)->c_str());
			if (cacheBufferMb == 0) {
				cacheBufferMb = 8;
			}
			params.erase(it, it + 2);
		}

		double filter = 1.0;
		uint32_t kmerLength = 18;

		it = find(params.begin(), std::prev(params.end()), "-f"); // minhash threshold
		if (it != std::prev(params.end())) {
			filter = atof(std::next(it)->c_str());
			params.erase(it, it + 2);
		}

		it = find(params.begin(), std::prev(params.end()), "-k"); // kmer length
		if (it != std::prev(params.end())) {
			kmerLength = atoi(std::next(it)->c_str());
			params.erase(it, it + 2);
		}

		// all modes need at least 2 parameters
		if (params.size() >= 2) {
			const string& mode = params[0];
			// building from kmers or genomes
			if (params.size() == 3 && mode == "build") {
				cout << "Database building mode (fasta genomes)" << endl;
				return runBuildDatabase(params[1], params[2], InputFile::GENOME, filter, kmerLength);
			} else 	if (params.size() == 3 && mode == "build-kmers") {
				cout << "Database building mode (kmers)" << endl;
				return runBuildDatabase(params[1], params[2], InputFile::KMC, filter, 0);
			} 
			else if (params.size() == 3 && (mode == "build-mh")) {
				cout << "Database building mode (mihashed kmers)" << endl;
				return runBuildDatabase(params[1], params[2], InputFile::MINHASH, 1.0, 0);
			}
			else if (params.size() == 3 && mode == "all2all") {
				cout << "All versus all comparison" << endl;
				return runAllVsAll(params[1], params[2]);
			}
			else if (params.size() == 4 && mode == "new2all") {
				cout << "Set of new samples (fasta genomes) versus entire database comparison" << endl;
				return runNewVsAll(params[1], params[2], params[3], InputFile::GENOME);
			}
			else if (params.size() == 4 && mode == "new2all-kmers") {
				cout << "Set of new samples (kmers) versus entire database comparison" << endl;
				return runNewVsAll(params[1], params[2], params[3], InputFile::KMC);
			}
			else if (params.size() == 4 && mode == "one2all") {
				cout << "One new sample (kmers) versus entire database comparison" << endl;
				return runOneVsAll(params[1], params[2], params[3]);
			}
			else if (params.size() == 3 && mode == "list-patterns") {
				cout << "Listing all patterns" << endl;
				return runListPatterns(params[1], params[2]);
			}
			else if (params.size() == 3 && mode == "minhash") {
				double filter = std::atof(params[1].c_str());
				if (filter > 0) {
					cout << "Minhashing k-mers" << endl;
					return runMinHash(params[2], filter);
				}
			}
			else if (params.size() == 2 && mode == "distance") {
				cout << "Calculating distance measures" << endl;
				return runDistanceCalculation(params[1]);
			}
			else if (params.size() == 3 && mode == "analyze") {
				cout << "Analyzing database" << endl;
				return runAnalyzeDatabase(params[1], params[2]);
			}
		}
	}

	showInstructions();
	return 0;
}

// *****************************************************************************************
//
int Console::runMinHash(const std::string& multipleKmcSamples, double filterValue) {
	cout << "Minhashing samples..." << endl;

	std::chrono::duration<double> loadingTime, processingTime;

	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<MinHashFilter>(filterValue, 0);

	Loader loader(filter, InputFile::KMC, numThreads);
	loader.configure(multipleKmcSamples);

	LOG_VERBOSE << "Init prefetch (group " << loader.getCurrentFileId() << ")" << endl;
	loader.initPrefetch();

	LOG_DEBUG << "Starting loop..." << endl;

	for (;;) {
		auto start = std::chrono::high_resolution_clock::now();
		LOG_VERBOSE << "Wait for prefetcher... " << endl;
		loader.waitForPrefetch();
		LOG_VERBOSE << "Init Loading... " << endl;
		loader.initLoad();
		LOG_VERBOSE << "Load... " << endl;
		loader.waitForLoad();
		loadingTime += std::chrono::high_resolution_clock::now() - start;

		start = std::chrono::high_resolution_clock::now();
		LOG_VERBOSE << "Init prefetch (group " << loader.getCurrentFileId() << ")" << endl;
		loader.initPrefetch();
		if (!loader.getLoadedTasks().size()) {
			break;
		}

		LOG_VERBOSE << "Process loaded tasks..." << endl;
		for (const auto& entry : loader.getLoadedTasks()) {
			auto task = entry.second;
			MihashedInputFile file(nullptr); 
			file.store(task->filePath, *task->kmers, task->kmerLength, filterValue);
		}

		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	return 0;
}

// *****************************************************************************************
//
int Console::runBuildDatabase(
	const std::string& multipleSamples, 
	const std::string dbFilename, 
	InputFile::Format inputFormat, 
	double filterValue,
	uint32_t kmerLength){

	cout << "Processing samples..." << endl;
	
	LOG_DEBUG << "Creating FastKmerDb object" << endl;
	FastKmerDb* db = new FastKmerDb(numThreads, cacheBufferMb);

	std::chrono::duration<double> loadingTime, processingTime;
	
	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<MinHashFilter>(filterValue, kmerLength);

	Loader loader(filter, inputFormat, numThreads);
	loader.configure(multipleSamples);

	loader.initPrefetch();

	LOG_DEBUG << "Starting loop..." << endl;
	auto totalStart = std::chrono::high_resolution_clock::now();
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
			db->addKmers(task->sampleName, *task->kmers, task->kmerLength, task->fraction);
			show_progress(*db);
		}
		
		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	loader.waitForPrefetch();

	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	cout << endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Processing time: " << processingTime.count() << endl
		<< db->printDetailedTimes() << endl
		<< "STATISTICS" << endl << db->printStats() << endl;

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

// *****************************************************************************************
//
int Console::runAllVsAll(const std::string& dbFilename, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	std::ofstream ofs(similarityFile);
	FastKmerDb* db = new FastKmerDb(numThreads, cacheBufferMb);;

	std::chrono::duration<double> dt;
	cout << "Loading k-mer database " << dbFilename << "...";
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db->deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl << db->printStats() << endl;


	cout << "Calculating matrix of common k-mers...";
	start = std::chrono::high_resolution_clock::now();
	LowerTriangularMatrix<uint32_t> matrix;
	db->calculateSimilarity(matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing matrix of common k-mers in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db->getKmerLength() << " fraction: " << db->getFraction() << " ,db-samples ,";
	std::copy(db->getSampleNames().cbegin(), db->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl << "query-samples,total-kmers,";
	std::copy(db->getSampleKmersCount().cbegin(), db->getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl;

	for (int sid = 0; sid < db->getSamplesCount(); ++sid) {
		ofs << db->getSampleNames()[sid] << ", " << db->getSampleKmersCount()[sid] << ", ";
		matrix.saveRow(sid, ofs);
		ofs << endl;
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Releasing memory...";
	start = std::chrono::high_resolution_clock::now();
	delete db;
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	return 0;
}

// *****************************************************************************************
//
int Console::runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads, cacheBufferMb);

	std::chrono::duration<double> dt;

	cout << "Loading k-mer database " << dbFilename << "...";
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl << db.printStats() << endl;
		
	cout << "Loading sample kmers...";

	start = std::chrono::high_resolution_clock::now();
	
	std::vector<kmer_t> queryKmers;
	std::vector<uint32_t> positions;
	uint32_t kmerLength;
	shared_ptr<MinHashFilter> filter = shared_ptr<MinHashFilter>(new MinHashFilter(db.getFraction(), db.getKmerLength()));
	
	KmcInputFile file(filter);
	double dummy;
	if (!file.open(singleKmcSample) || !file.load(queryKmers, positions, kmerLength, dummy)) {
		cout << "FAILED";
		return -1;
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl
		<< "Number of k-mers: " << queryKmers.size() << endl
		<< "Minhash fraction: " << db.getFraction() << endl;

	cout << "Calculating similarity vector...";
	start = std::chrono::high_resolution_clock::now();
	std::vector<uint32_t> sims;
	db.calculateSimilarity(queryKmers, sims);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing similarity vector in " << similarityFile << "...";
	std::ofstream ofs(similarityFile);

	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));

	ofs << endl << "query-samples,total-kmers,";
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl << singleKmcSample << "," << queryKmers.size() << ",";
	std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));

	ofs.close();
	cout << "OK" << endl;

	return 0;
}


// *****************************************************************************************
//
int Console::runNewVsAll(const std::string& dbFilename, const std::string& multipleSamples, const std::string& similarityFile, InputFile::Format inputFormat) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads, cacheBufferMb);
	
	std::chrono::duration<double> loadingTime, processingTime, dt;

	cout << "Loading k-mer database " << dbFilename << "...";
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl << db.printStats() << endl;

	// Opening file
	std::ofstream ofs(similarityFile);
	
	cout << "Storing matrix of common k-mers in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl << "query-samples,total-kmers,";
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl;

	cout << "Loading queries...";
	start = std::chrono::high_resolution_clock::now();

	LOG_DEBUG << "Creating Loader object..." << endl;
	shared_ptr<MinHashFilter> filter = shared_ptr<MinHashFilter>(new MinHashFilter(db.getFraction(), db.getKmerLength()));

	Loader loader(filter, inputFormat, numThreads);
	loader.configure(multipleSamples);

	loader.initPrefetch();

	LOG_DEBUG << "Starting loop..." << endl;

	std::vector<uint32_t> sims;

	auto totalStart = std::chrono::high_resolution_clock::now();
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

			// use position vector for storing common kmer counts
			sims.clear();
			db.calculateSimilarity(*task->kmers, sims);
			ofs << endl << task->sampleName << "," << task->kmers->size() << ",";
			std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));
		}

		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	loader.waitForPrefetch();

	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	cout << endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl
		<< "Loading k-mers: " << loadingTime.count() << endl
		<< "Processing time: " << processingTime.count() << endl;


	return 0;
}

// *****************************************************************************************
//
int Console::runDistanceCalculation(const std::string& similarityFilename) {

	std::vector<size_t> kmersCount;
	uint32_t kmerLength;
	
	cout << "Loading file with common k-mer counts: " << similarityFilename << "...";
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

	string tmp, in;
	double fraction;
	similarityFile >> tmp >> kmerLength >> tmp >> fraction >> tmp;

	getline(similarityFile, in); // copy sample names to output files
	for (auto & f : files) { f << "kmer-length: " << kmerLength << " fraction: " << fraction << in  << endl; }

	getline(similarityFile, in); // get number of kmers for all samples
	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss2(in);
	iss2 >> tmp >> tmp;
	std::copy(std::istream_iterator<size_t>(iss2), std::istream_iterator<size_t>(), std::back_inserter(kmersCount));

	for (int i = 0; getline(similarityFile, in); ++i) {
		std::replace(in.begin(), in.end(), ',', ' ');
		istringstream iss(in);
		uint64_t queryKmersCount = 0;
		string queryName;
		iss >> queryName >> queryKmersCount;
		
		for (int m = 0; m < metrics.size(); ++m) {
			files[m] << queryName << ",";
		}

		size_t intersection;
	
		for (int j = 0; iss >> intersection; ++j) {
			if (j >= kmersCount.size()) {
				cout << "Invalid file format!";
				return 0;
			}
			// calculate all metrices
			for (int m = 0; m < metrics.size(); ++m) {
				files[m] << metrics[m](intersection, queryKmersCount, kmersCount[j], kmerLength) << ","; 
			}
		}
		
		for (auto & f : files) { f << endl; }
	}

	cout << "OK" << endl;

	return 0;
}



// *****************************************************************************************
//
int Console::runAnalyzeDatabase(const std::string & multipleKmcSamples, const std::string & dbFilename)
{
	std::chrono::duration<double> loadingTime, processingTime;

	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads, cacheBufferMb);

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl << db.printStats();

	Analyzer aligner;
	
	aligner.printStats(db);

	std::shared_ptr<AbstractFilter> filter = aligner.selectSeedKmers(db, std::max((size_t)1, db.getSamplesCount() / 100));
	return 0;
	Loader loader(filter, InputFile::GENOME, numThreads, true);
	loader.configure(multipleKmcSamples);

	loader.initPrefetch();
	
	LOG_DEBUG << "Starting loop..." << endl;
	auto totalStart = std::chrono::high_resolution_clock::now();
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
			
			// use task result
		}

		loader.getLoadedTasks().clear();
		processingTime += std::chrono::high_resolution_clock::now() - start;
	}

	loader.waitForPrefetch();

	aligner(db);

	return 0;
}

// *****************************************************************************************
//
int Console::runListPatterns(const std::string& dbFilename, const std::string& patternFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	FastKmerDb db(numThreads, cacheBufferMb);

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl << db.printStats() << endl;
		
	cout << "Storing patterns in " << patternFile << "...";
	std::ofstream ofs(patternFile);
	db.savePatterns(ofs);
	ofs.close();
	cout << "OK" << endl; 

	return 0;
}

// *****************************************************************************************
//
void Console::showInstructions() {
	cout << "USAGE" << endl << endl

		<< "kmer-db <mode> [options] <positional arguments>" << endl << endl

		<< "Modes:" << endl
		<< "  build - building a database from genomes," << endl
		<< "  build-kmers - building a database from k-mers," << endl
		<< "  build-mh - building a database from minhashed k-mers," << endl
		<< "  minhash - minhashing k-mers," << endl
		<< "  all2all - calculating number of common k-mers between all samples in the database," << endl
		<< "  one2all - calculating number of common kmers between single sample and all the samples in the database," << endl
		<< "  distance - calculating similarities / distances." << endl
		<< "Global options:" << endl
		<< "  -t <threads> - number of threads (default: number of available cores)," << endl
		<< "  -buffer <size_mb> - size of cache buffer in megabytes, applies to all2all mode" << endl
		<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)." << endl
		<< "The meaning of the positional arguments depends on the selected mode." << endl << endl

		<< "Building k-mer database:" << endl
		<< "  kmer-db build [-f <filter> -k <kmer-length>] <sample_list> <database>" << endl
		<< "  kmer-db build-kmers [-f <filter>] <sample_list> <database>" << endl
		<< "  kmer_db build-mh <sample_list> <database>" << endl
		<< "    sample_list (input) - file containing list of samples, either raw KMC files (build mode) or minhashed (build-mh)" << endl
		<< "    database (output) - file with generated k-mer database" << endl
		<< "    -f <filter> - number from [0, 1] interval determining a fraction of all k-mers to be accepted by the minhash filter" << endl
		<< "    -k <kmer_length> - length of k-mers (default: 18)." << endl << endl

		<< "Minhashing k-mers:" << endl
		<< "  kmer_db minhash <sample_list> <fraction>" << endl
		<< "    sample_list (input) - file containing list of KMC k-mer files (raw)" << endl
		<< "    fraction (input) - fraction of kmers passing the filter" << endl << endl

		<< "Calculating number of common k-mers for all the samples in the database:" << endl
		<< "  kmer_db all2all <database> <common_matrix>" << endl
		<< "    database (input) - k-mer database file" << endl
		<< "    commony_matrix (output) - file containing matrix with numbers of common k-mers" << endl << endl

		<< "Calculating number of common kmers between single sample and all the samples in the database:" << endl
		<< "  kmer_db one2all <database> <sample> <common_vector>" << endl
		<< "    database (input) - k-mer database file" << endl
		<< "    sample (input) - sample name (corresponding KMC k-mer files must exist)" << endl
		<< "    common_vector (output) - file containing vector with numbers of common k-mers" << endl << endl

		<< "Calculating similarity/distance metrices (vectors) on the basis of common k-mers:" << endl
		<< "  kmer_db distance <common_file>" << endl
		<< "    common_file (input) - file containing matrix/vector with numbers of common k-mers" << endl
		<< "This mode generates set of files:" << endl
		<< "  <common_file>.jaccard" << endl
		<< "  <common_file>.min" << endl
		<< "  <common_file>.max" << endl
		<< "  <common_file>.cosine" << endl
		<< "  <common_file>.mash" << endl;
}

