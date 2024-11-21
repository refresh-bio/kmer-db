/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <algorithm>
#include <string>
#include <iostream>
#include <iterator>
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
#include "loader_ex.h"
#include "prefix_kmer_db.h"
#include "kmer_extract.h"

#include "../libs/refresh/sort/lib/pdqsort_par.h"

using namespace std;


void BuildConsole::run(const Params& params){
	if (params.files.size() != 2) {
		throw usage_error(params.mode);
	}

	LOG_NORMAL("Building database (from " << InputFile::format2string(params.inputFormat) << ")" << endl);
//	const std::string& multipleSamples(params.files[0]);
	const std::string multipleSamples(params.files[0]);
	const std::string dbFilename(params.files[1]);

	LOG_DEBUG("Creating PrefixKmerDb object" << endl);
	AbstractKmerDb* db = new PrefixKmerDb(params.numThreads);
	std::shared_ptr<MinHashFilter> filter;

	if (params.extendDb) {
		std::ifstream ifs;
		LOG_NORMAL("Loading k-mer database " << dbFilename << "...");
		ifs.open(dbFilename, std::ios::binary);
		if (!ifs || !db->deserialize(ifs)) {
			throw runtime_error("Cannot open k-mer database " + dbFilename);
		}
		filter = std::make_shared<MinHashFilter>(db->getFraction(), db->getStartFraction(), db->getKmerLength());
	}
	else {
		filter = std::make_shared<MinHashFilter>(params.fraction, params.fractionStart, params.kmerLength);
	}

	std::chrono::duration<double> sortingTime{ 0 }, processingTime{ 0 };
	
	LOG_NORMAL("Processing samples..." << endl);
	LOG_DEBUG("Creating Loader object..." << endl);

	LoaderEx loader(filter, params.inputFormat, params.numReaderThreads, params.numThreads, params.multisampleFasta);
	loader.configure(multipleSamples);

	LOG_DEBUG("Starting loop..." << endl);
	auto totalStart = std::chrono::high_resolution_clock::now();
	int sample_id = 0;
	for (; !loader.isCompleted(); ++sample_id) {
		auto partialTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);
		LOG_VERBOSE("Processing time: " << partialTime.count() << ", loader buffers: " << (loader.getBytes() >> 20) << " MB" << endl);
		
		auto task = loader.popTask(sample_id);

		if (task) {
			if ((sample_id + 1) % 10 == 0) {
				size_t cnt = loader.getSamplesCount();
				if (cnt > 0) {
					LOG_NORMAL("\r" << sample_id + 1 << "/" << cnt << "...");
				}
				else {
					LOG_NORMAL("\r" << sample_id + 1 << "...");
				}
			}

			auto start = std::chrono::high_resolution_clock::now();

			// postprocess k-mers if neccessary
			if (params.inputFormat == InputFile::Format::GENOME) {
//				KmerHelper::sortAndUnique(task->kmers, task->kmersCount, params.numThreads);

				refresh::sort::pdqsort_branchless_tp(refresh::sort::pdqsort_adjust_threads(task->kmersCount, params.numThreads), task->kmers, task->kmers + task->kmersCount, atp);
//				refresh::sort::pdqsort_branchless(refresh::sort::pdqsort_adjust_threads(task->kmersCount, params.numThreads), task->kmers, task->kmers + task->kmersCount);
//				refresh::sort::pdqsort_branchless(task->kmers, task->kmers + task->kmersCount);
//				stable_sort(task->kmers, task->kmers + task->kmersCount);
				auto it = std::unique(task->kmers, task->kmers + task->kmersCount);
				task->kmersCount = it - task->kmers;
			}
			else if (params.inputFormat == InputFile::Format::KMC) {
				KmerHelper::sort(task->kmers, task->kmersCount, params.numThreads);
			}

			auto start2 = std::chrono::high_resolution_clock::now();
			sortingTime += start2 - start;

			db->addKmers(task->sampleName, task->kmers, task->kmersCount, task->kmerLength, task->fraction, atp);
			processingTime += std::chrono::high_resolution_clock::now() - start2;
			
			loader.releaseTask(*task);
			LOG_VERBOSE(db->printProgress() << endl);
		}
	}

	LOG_NORMAL("\r" << sample_id << "/" << sample_id << "                      " << endl);

	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	LOG_NORMAL(endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl
		<< "Kmer sorting/unique time: " << sortingTime.count() << endl
		<< "Database update time:" << processingTime.count() << endl);
#ifdef COLLECT_DETAILED_TIMES
	LOG_NORMAL(db->printDetailedTimes() << endl);
#endif
	LOG_NORMAL("STATISTICS" << endl << db->printStats() << endl);

	std::chrono::duration<double> dt{ 0 };

	LOG_NORMAL("Serializing database..." << endl);

	const size_t io_buffer_size = 64 << 20;
	std::ofstream ofs;
	ofs.open(dbFilename, ios_base::out | std::ios::binary);
	char* io_buffer = new char[io_buffer_size];
	ofs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);
	db->serialize(ofs, true);
	ofs.close();
	delete[] io_buffer;

	LOG_NORMAL(endl << "Releasing memory...");
	auto start = std::chrono::high_resolution_clock::now();
	delete db;
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL("OK (" << dt.count() << " seconds)" << endl);
}