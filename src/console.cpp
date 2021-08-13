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

#ifndef __APPLE__
#include <omp.h>
#endif

#include "console.h"
#include "kmer_db.h"
#include "loader_ex.h"
#include "version.h"
#include "analyzer.h"
#include "similarity_calculator.h"
#include "prefix_kmer_db.h"
#include "kmer_extract.h"


using namespace std;



const string Params::MODE_BUILD = "build";
const string Params::MODE_MINHASH = "minhash";
const string Params::MODE_ALL_2_ALL = "all2all";
const string Params::MODE_NEW_2_ALL = "new2all";
const string Params::MODE_ONE_2_ALL = "one2all";
const string Params::MODE_DISTANCE = "distance";

const string Params::SWITCH_KMC_SAMPLES = "-from-kmers";
const string Params::SWITCH_MINHASH_SAMPLES = "-from-minhash";
const string Params::SWITCH_MULTISAMPLE_FASTA = "-multisample-fasta";
const string Params::SWITCH_PHYLIP_OUT = "-phylip-out";
const string Params::SWITCH_EXTEND_DB = "-extend";
const string Params::SWITCH_SPARSE = "-sparse";

const string Params::OPTION_FRACTION = "-f";
const string Params::OPTION_FRACTION_START = "-f-start";
const string Params::OPTION_LENGTH = "-k";
const string Params::OPTION_VERBOSE = "-v";
const string Params::OPTION_DEBUG = "-vv";
const string Params::OPTION_THREADS = "-t";
const string Params::OPTION_READER_THREADS = "-rt";
const string Params::OPTION_BUFFER = "-buffer";


// *****************************************************************************************
//
Console::Console() {
	availableMetrics["jaccard"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / (cnt1 + cnt2 - common); };
	availableMetrics["min"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::min(cnt1, cnt2); };
	availableMetrics["max"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::max(cnt1, cnt2); };
	availableMetrics["cosine"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / sqrt(cnt1 * cnt2); };
	availableMetrics["mash"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / (queryCnt + dbCnt - common);
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1)); 
	};

	availableMetrics["mash-query"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / queryCnt;
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
	};
}


// *****************************************************************************************
//
int Console::parse(int argc, char** argv) {

	cout << "Kmer-db version " << VERSION << " (" << DATE << ")" << endl 
		<< "S. Deorowicz, A. Gudys, M. Dlugosz, M. Kokot, and A. Danek (c) 2018" << endl << endl;
	
	time_t rawtime;
	struct tm * timeinfo;

	// set initial values of parameters
	numThreads = 0;
	numReaderThreads = 0;
	cacheBufferMb = 8;
	multisampleFasta = false;
	fraction = 1.0;
	fractionStart = 0.0;
	kmerLength = 18;
	sparse = false;

	InputFile::Format inputFormat = InputFile::GENOME;
	bool extendDb = false;
	
	std::vector<string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* c)->string { return c; });

	int status = -1;

	if (params.size()) {

		// search for switches and options
		if (findSwitch(params, Params::OPTION_VERBOSE)) { // verbose mode
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
		}

		if (findSwitch(params, Params::OPTION_DEBUG)) { // verbose mode
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
			Log::getInstance(Log::LEVEL_DEBUG).enable();
		}

		multisampleFasta = findSwitch(params, Params::SWITCH_MULTISAMPLE_FASTA);
	
		findOption(params, Params::OPTION_FRACTION, fraction);				// minhash fraction
		findOption(params, Params::OPTION_FRACTION_START, fractionStart);	// minhash fraction start value
		findOption(params, Params::OPTION_LENGTH, kmerLength);				// kmer length
		
		findOption(params, Params::OPTION_THREADS, numThreads);			// number of threads
		if (numThreads <= 0) {
			numThreads = std::thread::hardware_concurrency();
		}

		// limit number of threads used in parallel sorts
#ifndef __APPLE__
		omp_set_num_threads(numThreads);
#endif
		
		findOption(params, Params::OPTION_READER_THREADS, numReaderThreads);	// number of threads
		if (numReaderThreads <= 0) {
			// more reader threads for smaller filters (from t/8 up to t)
			int invFraction = (int)(1.0 / fraction);
			numReaderThreads = std::max(std::min(numThreads, (numThreads / 8) * invFraction), 1); 
		}

		findOption(params, Params::OPTION_BUFFER, cacheBufferMb);	// size of temporary buffer in megabytes
		if (cacheBufferMb <= 0) {
			cacheBufferMb = 8;
		}

		if (findSwitch(params, Params::SWITCH_KMC_SAMPLES)) {
			inputFormat = InputFile::KMC;
			kmerLength = 0;
		}

		if (findSwitch(params, Params::SWITCH_MINHASH_SAMPLES)) {
			if (inputFormat == InputFile::KMC) {
				cout << "Error: " << Params::SWITCH_KMC_SAMPLES << " and " << Params::SWITCH_MINHASH_SAMPLES << " switches exclude one another." << endl;
				return 0;
			}
			inputFormat = InputFile::MINHASH;
			fraction = 1.0;
			kmerLength = 0;
		}

		sparse = findSwitch(params, Params::SWITCH_SPARSE);
		extendDb = findSwitch(params, Params::SWITCH_EXTEND_DB);
		
		// all modes need at least 2 parameters
		if (params.size() >= 2) {
			const string& mode = params[0];

			// detect obsolete modes
			if (mode == "build-kmers" || mode == "build-mh") {
				cout << "Error: build-kmers/build-mh modes are obsolete, use " << Params::SWITCH_KMC_SAMPLES << "/" << Params::SWITCH_MINHASH_SAMPLES << " switches instead." << endl;
				return 0;
			}

			time(&rawtime);
			timeinfo = localtime(&rawtime);
			cout << "Analysis started at " << asctime(timeinfo) << endl;

			// main modes
			if (params.size() == 3 && mode == Params::MODE_BUILD) {
				cout << "Database building mode (from " << InputFile::format2string(inputFormat) << ")" << endl;
				status = runBuildDatabase(params[1], params[2], inputFormat, extendDb);
			}
			else if (params.size() == 3 && mode == Params::MODE_ALL_2_ALL) {
				cout << "All versus all comparison" << endl;
				status = runAllVsAll(params[1], params[2]);
			}
			else if (params.size() == 4 && mode == Params::MODE_NEW_2_ALL) {
				cout << "Set of new samples  (from " << InputFile::format2string(inputFormat) << ") versus entire database comparison" << endl;
				
				status = runNewVsAll(params[1], params[2], params[3], inputFormat);
			}
			else if (params.size() == 4 && mode == Params::MODE_ONE_2_ALL) {
				cout << "One new sample  (from " << InputFile::format2string(inputFormat) << ") versus entire database comparison" << endl;
				status = runOneVsAll(params[1], params[2], params[3], inputFormat);
			}
			else if (params.size() == 3 && mode == Params::MODE_MINHASH) {
				fraction = std::atof(params[1].c_str());
				if (fraction > 0) {
					cout << "Minhashing k-mers" << endl;
					status = runMinHash(params[2], inputFormat);
				}
			}
			else if (params.size() >= 2 && mode == Params::MODE_DISTANCE) {

				bool phylipOut = findSwitch(params, Params::SWITCH_PHYLIP_OUT);

				// check selected metrics
				std::vector<string> metricNames;
				for (const auto& entry : availableMetrics) {
					if (findSwitch(params, entry.first)) {
						metricNames.push_back(entry.first);
					}
				}

				// if empty, add jacard
				if (metricNames.empty()) {
					metricNames.push_back("jaccard");
				}
				
				cout << "Calculating distance measures" << endl;
				status = runDistanceCalculation(params[1], metricNames, phylipOut);
			}
			// debug modes
			else if (params.size() == 3 && mode == "list-patterns") {
				cout << "Listing all patterns" << endl;
				status = runListPatterns(params[1], params[2]);
			}
			else if (params.size() == 3 && mode == "analyze") {
				cout << "Analyzing database" << endl;
				status = runAnalyzeDatabase(params[1], params[2]);
			}
		}
	}

	if (status == -1) {
		showInstructions();
	}
	else {
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cout << endl << "Analysis finished at " << asctime(timeinfo) << endl;
	}

	
	return 0;
}

// *****************************************************************************************
//
int Console::runMinHash(const std::string& multipleKmcSamples, InputFile::Format inputFormat) {
	cout << "Minhashing samples..." << endl;

	std::chrono::duration<double> loadingTime{ 0 }, processingTime{ 0 };

	LOG_DEBUG << "Creating Loader object..." << endl;

	auto filter = std::make_shared<MinHashFilter>(fraction, 0, kmerLength);

	LoaderEx loader(filter, inputFormat, numReaderThreads, numThreads, multisampleFasta);
	loader.configure(multipleKmcSamples);

	LOG_DEBUG << "Starting loop..." << endl;
	auto totalStart = std::chrono::high_resolution_clock::now();
	for (int i = 0; !loader.isCompleted(); ++i) {
		auto partialTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);
		LOG_VERBOSE << "Processing time: " << partialTime.count() << ", loader buffers: " << (loader.getBytes() >> 20) << " MB" << endl;

		auto task = loader.popTask(i);

		if (task) {
			auto start = std::chrono::high_resolution_clock::now();

			MihashedInputFile file(nullptr);

			// postprocess k-mers if neccessary
			if (inputFormat == InputFile::Format::GENOME) {
				KmerHelper::sortAndUnique(task->kmers, task->kmersCount);
			}
			else if (inputFormat == InputFile::Format::KMC) {
				KmerHelper::sort(task->kmers, task->kmersCount);
			}

			file.store(task->filePath, task->kmers, task->kmersCount, task->kmerLength, fraction);

			processingTime += std::chrono::high_resolution_clock::now() - start;
			loader.releaseTask(*task);
		}
	}

	return 0;
}

// *****************************************************************************************
//
int Console::runBuildDatabase(
	const std::string& multipleSamples, 
	const std::string dbFilename, 
	InputFile::Format inputFormat,
	bool extendDb){

	LOG_DEBUG << "Creating PrefixKmerDb object" << endl;
	AbstractKmerDb* db = new PrefixKmerDb(numThreads);
	std::shared_ptr<MinHashFilter> filter;

	if (extendDb) {
		std::ifstream ifs;
		cout << "Loading k-mer database " << dbFilename << "..." ;
		ifs.open(dbFilename, std::ios::binary);
		if (!ifs || !db->deserialize(ifs)) {
			cout << "FAILED!" << endl;
			return -1;
		}
		filter = std::make_shared<MinHashFilter>(db->getFraction(), db->getStartFraction(), db->getKmerLength());
	}
	else {
		filter = std::make_shared<MinHashFilter>(fraction, fractionStart, kmerLength);
	}

	std::chrono::duration<double> sortingTime{ 0 }, processingTime{ 0 };
	
	cout << "Processing samples..." << endl;
	LOG_DEBUG << "Creating Loader object..." << endl;

	LoaderEx loader(filter, inputFormat, numReaderThreads, numThreads, multisampleFasta);
	loader.configure(multipleSamples);

	LOG_DEBUG << "Starting loop..." << endl;
	auto totalStart = std::chrono::high_resolution_clock::now();
	int sample_id = 0;
	for (; !loader.isCompleted(); ++sample_id) {
		auto partialTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);
		LOG_VERBOSE << "Processing time: " << partialTime.count() <<  ", loader buffers: " << (loader.getBytes() >> 20) << " MB" << endl;
		
		auto task = loader.popTask(sample_id);

		if (task) {
			if ((sample_id + 1) % 10 == 0) {
				size_t cnt = loader.getSamplesCount();
				if (cnt > 0) {
					cout << "\r" << sample_id + 1 << "/" << cnt << "..." << std::flush;
				}
				else {
					cout << "\r" << sample_id + 1 << "..." << std::flush;
				}
			}

			auto start = std::chrono::high_resolution_clock::now();

			// postprocess k-mers if neccessary
			if (inputFormat == InputFile::Format::GENOME) {
				KmerHelper::sortAndUnique(task->kmers, task->kmersCount);
			}
			else if (inputFormat == InputFile::Format::KMC) {
				KmerHelper::sort(task->kmers, task->kmersCount);
			}
			sortingTime += std::chrono::high_resolution_clock::now() - start;

			start = std::chrono::high_resolution_clock::now();
			db->addKmers(task->sampleName, task->kmers, task->kmersCount, task->kmerLength, task->fraction);
			processingTime += std::chrono::high_resolution_clock::now() - start;
			
			loader.releaseTask(*task);
			LOG_VERBOSE << db->printProgress() << endl;
		}
	}

	cout << "\r" << sample_id << "/" << sample_id << "                      " << endl << std::flush;

	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	cout << endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl
		<< "Kmer sorting/unique time: " << sortingTime.count() << endl
		<< "Database update time:" << processingTime.count() << endl
		<< db->printDetailedTimes() << endl
		<< "STATISTICS" << endl << db->printStats() << endl;

	std::chrono::duration<double> dt{ 0 };

	cout << "Serializing database..." << endl;
	std::ofstream ofs;
	ofs.open(dbFilename, std::ios::binary);
	db->serialize(ofs, true);
	ofs.close();

	cout << endl << "Releasing memory...";
	auto start = std::chrono::high_resolution_clock::now();
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
	PrefixKmerDb* db = new PrefixKmerDb(numThreads);
	SimilarityCalculator calculator(numThreads, cacheBufferMb);
	
	std::chrono::duration<double> dt{ 0 };
	cout << "Loading k-mer database " << dbFilename << "..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db->deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	
	cout << "Calculating matrix of common k-mers...";
	start = std::chrono::high_resolution_clock::now();
	LowerTriangularMatrix<uint32_t> matrix;
	calculator.all2all(*db, matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing matrix of common k-mers in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db->getKmerLength() << " fraction: " << db->getFraction() << " ,db-samples ,";
	std::copy(db->getSampleNames().cbegin(), db->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl << "query-samples,total-kmers,";
	std::copy(db->getSampleKmersCount().cbegin(), db->getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl;

	for (size_t sid = 0; sid < db->getSamplesCount(); ++sid) {
		ofs << db->getSampleNames()[sid] << ", " << db->getSampleKmersCount()[sid] << ", ";
		if (sparse) {
			matrix.saveRowSparse(sid, ofs);
		}
		else {
			matrix.saveRow(sid, ofs);
		}
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
int Console::runOneVsAll(const std::string& dbFilename, const std::string& sampleFasta, const std::string& similarityFile, InputFile::Format inputFormat) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(numThreads);
	SimilarityCalculator calculator(numThreads, cacheBufferMb);

	std::chrono::duration<double> dt{ 0 };

	cout << "Loading k-mer database " << dbFilename << ":" << endl;
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl << db.printStats() << endl;
		
	cout << "Loading sample kmers...";

	start = std::chrono::high_resolution_clock::now();
	
	std::vector<kmer_t> kmersBuffer;
	std::vector<uint32_t> positions;
	uint32_t kmerLength;
	shared_ptr<MinHashFilter> filter = shared_ptr<MinHashFilter>(new MinHashFilter(db.getFraction(), db.getStartFraction(), db.getKmerLength()));
	
	std::shared_ptr<InputFile> file;
	
	if (inputFormat == InputFile::KMC) {
		file = std::make_shared<KmcInputFile>(filter->clone());
	}
	else if (inputFormat == InputFile::MINHASH) {
		file = std::make_shared<MihashedInputFile>(filter->clone());
	}
	else {
		file = std::make_shared<GenomeInputFile>(filter->clone(), false);
	}

	double dummy;
	size_t queryKmersCount;
	kmer_t* queryKmers;

	if (!file->open(sampleFasta) || !file->load(kmersBuffer, positions, queryKmers, queryKmersCount, kmerLength, dummy)) {
		cout << "FAILED";
		return -1;
	}

	// postprocess k-mers if neccessary
	if (inputFormat == InputFile::Format::GENOME) {
		KmerHelper::sortAndUnique(queryKmers, queryKmersCount);
	}

	if (kmerLength != db.getKmerLength()) {
		cout << "Error: sample and database k-mer length differ." << endl;
		return -1;
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl
		<< "Number of k-mers: " << queryKmersCount << endl
		<< "Minhash fraction: " << db.getFraction() << endl;

	cout << "Calculating similarity vector...";
	start = std::chrono::high_resolution_clock::now();
	std::vector<uint32_t> sims;
	calculator.one2all(db, queryKmers, queryKmersCount, sims);
	dt = std::chrono::high_resolution_clock::now() - start;
	cout << "OK (" << dt.count() << " seconds)" << endl;

	cout << "Storing similarity vector in " << similarityFile << "...";
	std::ofstream ofs(similarityFile);

	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));

	ofs << endl << "query-samples,total-kmers,";
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl << sampleFasta << "," << queryKmersCount << ",";
	std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));

	ofs.close();
	cout << "OK" << endl;

	return 0;
}


// *****************************************************************************************
//
int Console::runNewVsAll(const std::string& dbFilename, const std::string& multipleSamples, const std::string& similarityFile, InputFile::Format inputFormat) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(numThreads);
	SimilarityCalculator calculator(numThreads, cacheBufferMb);

	std::chrono::duration<double> loadingTime{ 0 }, processingTime{ 0 }, dt{ 0 };

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
	
	cout << "Storing matrix of common k-mers in " << similarityFile << endl;
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl << "query-samples,total-kmers,";
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	
	cout << "Loading queries...";
	start = std::chrono::high_resolution_clock::now();

	LOG_DEBUG << "Creating Loader object..." << endl;
	shared_ptr<MinHashFilter> filter = shared_ptr<MinHashFilter>(new MinHashFilter(db.getFraction(), db.getStartFraction(), db.getKmerLength()));

	LoaderEx loader(filter, inputFormat, numReaderThreads, numThreads, multisampleFasta);
	loader.configure(multipleSamples);
	cout << endl;

	std::vector<uint32_t> sims;
	
	cout << "Processing queries..." << endl;
	auto totalStart = std::chrono::high_resolution_clock::now();

	// create set of buffers for storing similarities
	std::vector<std::vector<uint32_t>> buffers(loader.getOutputBuffersCount());
	RegisteringQueue<int> freeBuffersQueue(1);
	for (size_t i = 0; i < buffers.size(); ++i) {
		freeBuffersQueue.Push(i);
	}

	SynchronizedPriorityQueue<std::shared_ptr<SampleTask>> similarityQueue(numThreads);
	std::vector<thread> workers(numThreads);
	std::atomic<int> sample_id{ 0 };

	for (int tid = 0; tid < numThreads; ++tid) {
		workers[tid] = thread([&db, &loader, &freeBuffersQueue, &similarityQueue, &buffers, &calculator, &sample_id, tid]() {
			while (!loader.isCompleted()) {
				int task_id = sample_id.fetch_add(1);
				auto task = loader.popTask(task_id);
				
				if (task) {
					task->bufferId2 = -1;
					freeBuffersQueue.Pop(task->bufferId2);
					buffers[task->bufferId2].clear();
					
					// only unique k-mers are needed
					KmerHelper::unique(task->kmers, task->kmersCount);
					
					calculator.one2all<false>(db, task->kmers, task->kmersCount, buffers[task->bufferId2]);
					similarityQueue.Push(task_id, task);
				
					LOG_DEBUG << "(" << task_id + 1 << ") -> similarity queue, tid:"  << tid << ", buf:" << task->bufferId2 << endl;
				}
			}

			similarityQueue.MarkCompleted();
			LOG_DEBUG << "processing finished, tid: " << tid << endl;
		});
	}

	// gather results in one thread
	for (int task_id = 0; !similarityQueue.IsCompleted(); ++task_id) {
		
		if ((task_id + 1) % 10 == 0) {
			cout << "\r" << task_id + 1 << "...                      " << std::flush;
		}
		
		std::shared_ptr<SampleTask> task;
		similarityQueue.Pop(task_id, task);
		LOG_DEBUG << "similarity queue -> (" << task_id + 1 << ", " << task->sampleName << "), buf:" << task->bufferId2 << endl;
		const auto& row = buffers[task->bufferId2];
		ofs << endl << task->sampleName << "," << task->kmersCount << ",";
		if (sparse) {
			for (int i = 0; i < row.size(); ++i) {
				if (row[i] > 0) {
					ofs << (i + 1) << ":" << row[i] << ",";
				}
			}
		}
		else {
			std::copy(row.begin(), row.end(), ostream_iterator<uint32_t>(ofs, ","));
		}
		freeBuffersQueue.Push(task->bufferId2);
		loader.releaseTask(*task);
	}

	// make sure all threads have finished
	for (auto &w : workers) {
		w.join();
	}

 	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	cout << endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl;


	return 0;
}

// *****************************************************************************************
//
int Console::runDistanceCalculation(const std::string& similarityFilename, const std::vector<string>& metricNames, bool usePhylip) {

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

	std::vector<metric_fun_t> metrics; 
	for (const auto& name : metricNames) {
		metrics.push_back(availableMetrics[name]);
	}
	
	std::vector<std::ofstream> files(metricNames.size());
//	std::vector<std::ofstream> histoFiles(metricNames.size());
	
	for (size_t i = 0; i < files.size(); ++i) {
		files[i].open(similarityFilename + "." + metricNames[i]);
//		histoFiles[i].open(similarityFilename + "." + metricNames[i] + ".histo");
	}

	std::vector<std::vector<size_t>> histograms(metrics.size(), std::vector<size_t>(100));

	string tmp, in;
	double fraction;
	similarityFile >> tmp >> kmerLength >> tmp >> fraction >> tmp;

	getline(similarityFile, in); // copy sample names to output files
	
	if (!usePhylip) {
		for (auto & f : files) {
			f << "kmer-length: " << kmerLength << " fraction: " << fraction << in << endl;
		}
	}

	getline(similarityFile, in); // get number of kmers for all samples
	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss2(in);
	iss2 >> tmp >> tmp;
	std::copy(std::istream_iterator<size_t>(iss2), std::istream_iterator<size_t>(), std::back_inserter(kmersCount));

	if (usePhylip) {
		for (auto & f : files) {
			f << kmersCount.size() << endl;
		}
	}

	std::vector<size_t> intersections(kmersCount.size());
	std::vector<double> values(kmersCount.size());

	char* outBuffer = new char[10000000];

	cout << "Processing rows..." << endl;
	for (int i = 0; getline(similarityFile, in); ++i) {
		if ((i + 1) % 10 == 0) {
			cout << "\r" << i + 1 << "/" << kmersCount.size() << "...";
		}
		// find first comma 
		auto pos = in.find(',');
		string queryName(in.begin(), in.begin() + pos);

		std::replace(in.begin(), in.end(), ',', ' ');
		istringstream iss(in);
		iss.seekg(pos + 1); // move right after first comma separator
		uint64_t queryKmersCount = 0;

		iss >> queryKmersCount;
	
		auto newEnd = std::copy(std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>(), intersections.begin());
		size_t numVals = newEnd - intersections.begin();

		for (size_t m = 0; m < metrics.size(); ++m) {
			auto& metric = metrics[m];
			
			std::transform(intersections.begin(), newEnd, kmersCount.begin(), values.begin(),
				[&metric, queryKmersCount, kmerLength](size_t intersection, size_t dbKmerCount)->double { return  metric(intersection, queryKmersCount, dbKmerCount, kmerLength); });
			
			char* ptr = outBuffer;
			memcpy(ptr, queryName.c_str(), queryName.size());
			ptr += queryName.size();
			*ptr = usePhylip ? ' ' : ',';
			++ptr;

			for (size_t j = 0; j < numVals; ++j) {
				
				if (values[j] == 0) {
					*ptr = '0';
					++ptr;
				}
				else {
					ptr += NumericConversions::Double2PChar(values[j], 6, ptr);
				}
					
//				int binId = (size_t)(values[j] * 100);
//				++histograms[m][binId];
				*ptr = usePhylip ? ' ' : ',';
				++ptr;
			}
			*ptr = 0;
			size_t len = ptr - outBuffer;
			files[m].write(outBuffer, len);
			files[m] << endl;
		}
	}

	for (size_t m = 0; m < metrics.size(); ++m) {
//		std::copy(histograms[m].begin(), histograms[m].end(), std::ostream_iterator<size_t>(histoFiles[m], ","));
	}

	cout << "OK" << endl;

	delete[] outBuffer;

	return 0;
}


// *****************************************************************************************
//
int Console::runAnalyzeDatabase(const std::string & multipleKmcSamples, const std::string & dbFilename)
{
	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(numThreads);

	cout << "Loading k-mer database " << dbFilename << "...";
	if (!dbFile || !db.deserialize(dbFile)) {
		cout << "FAILED";
		return -1;
	}
	cout << "OK" << endl << db.printStats();

	Analyzer analyzer;
	
	analyzer.printStats(db);

	std::shared_ptr<AbstractFilter> filter = analyzer.selectSeedKmers(db, std::max((size_t)1, db.getSamplesCount() / 100));
	return 0;
	LoaderEx loader(filter, InputFile::GENOME, numReaderThreads, numThreads, true);
	int numSamples = loader.configure(multipleKmcSamples);

	LOG_DEBUG << "Starting loop..." << endl;
	for (int i = 0; i < numSamples; ++i) {
		
		auto task = loader.popTask(i);
		// use task result
	}

	analyzer(db);

	return 0;
}

// *****************************************************************************************
//
int Console::runListPatterns(const std::string& dbFilename, const std::string& patternFile) {
	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(numThreads);

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
		<< "  " << Params::MODE_BUILD << " - building a database from genomes, k-mers, or minhashed k-mers," << endl

		<< "  " << Params::MODE_ALL_2_ALL << " - counting common k-mers - all samples in the database," << endl
		<< "  " << Params::MODE_NEW_2_ALL << " - counting common k-mers - set of new samples versus database," << endl
		<< "  " << Params::MODE_ONE_2_ALL << " - counting common k-mers - single sample versus database," << endl
		<< "  " << Params::MODE_DISTANCE << " - calculating similarities/distances." << endl
		<< "  " << Params::MODE_MINHASH << " - storing minhashed k-mers," << endl
		<< "Common options:" << endl
		<< "  " << Params::OPTION_THREADS << " <threads> - number of threads (default: number of available cores)," << endl
		<< "The meaning of other options and positional arguments depends on the selected mode." << endl << endl

		<< "Building a database:" << endl
		<< "  kmer-db " << Params::MODE_BUILD << " [" << Params::OPTION_LENGTH << " <kmer-length>] [" << Params::OPTION_FRACTION << " <fraction>] [" 
			<< Params::SWITCH_MULTISAMPLE_FASTA << "] [" << Params::SWITCH_EXTEND_DB << "] <sample_list> <database>" << endl
		<< "  kmer-db " << Params::MODE_BUILD << " " << Params::SWITCH_KMC_SAMPLES << " [" << Params::OPTION_FRACTION << " <fraction>] [" << Params::SWITCH_EXTEND_DB << "] <sample_list> <database>" << endl
		<< "  kmer-db " << Params::MODE_BUILD << " " << Params::SWITCH_MINHASH_SAMPLES << " [" << Params::SWITCH_EXTEND_DB << "] <sample_list> <database>" << endl
		<< "    sample_list (input) - file containing list of samples in one of the following formats:" << endl
		<< "                          fasta genomes or reads (default), KMC k-mers (" << Params::SWITCH_KMC_SAMPLES << "), or minhashed k-mers (" << Params::SWITCH_MINHASH_SAMPLES << ")," << endl
		<< "    database (output) - file with generated k-mer database," << endl
		<< "    " << Params::OPTION_LENGTH << " <kmer_length> - length of k-mers (default: 18)," << endl
		<< "    " << Params::OPTION_FRACTION << " <fraction> - fraction of all k-mers to be accepted by the minhash filter (default: 1)," << endl
		<< "    " << Params::SWITCH_MULTISAMPLE_FASTA << " - each sequence in a genome FASTA file is treated as a separate sample," << endl
		<< "	" << Params::SWITCH_EXTEND_DB << " - extend the existing database with new samples." << endl <<endl

		<< "Counting common k-mers for all the samples in the database:" << endl
		<< "  kmer-db " << Params::MODE_ALL_2_ALL << " [" << Params::OPTION_BUFFER << " <size_mb>] <database> <common_table>" << endl
		<< "    database (input) - k-mer database file." << endl
		<< "    common_table (output) - comma-separated table with number of common k-mers." << endl
		<< "    " << Params::OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes, applies to all2all mode" << endl
		<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)." << endl << endl

		<< "Counting common kmers between set of new samples and all the samples in the database:" << endl
		<< "  kmer-db " << Params::MODE_NEW_2_ALL << " [" << Params::SWITCH_MULTISAMPLE_FASTA << " | " << Params::SWITCH_KMC_SAMPLES << " | " << Params::SWITCH_MINHASH_SAMPLES << "] <database> <sample_list> <common_table>" << endl
		<< "    database (input) - k-mer database file." << endl
		<< "    sample_list (input) - file containing list of query samples in one of the supported formats (see build mode)." << endl
		<< "    common_table (output) - comma-separated table with number of common k-mers." << endl << endl

		<< "Counting common kmers between single sample and all the samples in the database:" << endl
		<< "  kmer-db " << Params::MODE_ONE_2_ALL << " [" << Params::SWITCH_MULTISAMPLE_FASTA << " | " << Params::SWITCH_KMC_SAMPLES << " | " << Params::SWITCH_MINHASH_SAMPLES << "] <database> <sample> <common_table>" << endl
		<< "    database (input) - k-mer database file." << endl
		<< "    sample (input) - query sample in one of the supported formats (see build mode)." << endl
		<< "    common_table (output) - comma-separated table with number of common k-mers." << endl << endl

		<< "Calculating similarities/distances on the basis of common k-mers:" << endl
		<< "  kmer-db " << Params::MODE_DISTANCE << " [<measures>] [" << Params::SWITCH_PHYLIP_OUT << "] <common_table>" << endl
		<< "    common_table (input) - comma-separated table with number of common k-mers." << endl
		<< "    " << Params::SWITCH_PHYLIP_OUT << " - store output distance matrix in a Phylip format" << endl	
		<< "    measures - names of the similarity/distance measures to be calculated, one or more of the following" << endl
		<< "               jaccard (default), min, max, cosine, mash." << endl

		<< "This mode generates a file with similarity/distance table for each selected measure." << endl
		<< "Name of the output file is produced by adding to the input file an extension with a measure name." << endl << endl

		<< "Storing minhashed k-mers:" << endl
		<< "  kmer-db " << Params::MODE_MINHASH << " [" << Params::OPTION_LENGTH << " <kmer-length>]" << " [" << Params::SWITCH_MULTISAMPLE_FASTA << "] <fraction> <sample_list>" << endl
		<< "  kmer-db " << Params::MODE_MINHASH << " " << Params::SWITCH_KMC_SAMPLES << " <fraction> <sample_list>" << endl
		<< "    fraction (input) - fraction of all k-mers to be accepted by the minhash filter." << endl
		<< "    sample_list (input) - file containing list of query samples in one of the supported formats (see build mode)." << endl
		<< "For each sample from the list, a binary file with *.minhash* extension containing filtered k-mers is created." << endl << endl;
		
}

