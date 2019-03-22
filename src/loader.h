#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "kmc_file_wrapper.h"
#include "queue.h"
#include "filter.h"

#include <vector>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>

// *****************************************************************************************
//
struct Task {
	size_t fileId;
	size_t threadId;
	std::string filePath;
	std::string sampleName;
	std::shared_ptr<InputFile> file;
	std::vector<kmer_t>* kmers;
	std::vector<uint32_t>* positions;
	uint32_t kmerLength;
	double fraction;

	// *****************************************************************************************
	//
	Task(size_t fileId, size_t threadId, const std::string& filePath) :
		fileId(fileId), threadId(threadId), filePath(filePath), file(nullptr), kmers(nullptr), positions(nullptr) {
	
		size_t pos = filePath.find_last_of("/\\");
		if (pos != string::npos) {
			sampleName = filePath.substr(pos + 1);
		}
		else {
			sampleName = filePath;
		}
	}
};

// *****************************************************************************************
//
class Loader {
public:
	
	size_t getCurrentFileId() const { return currentFileId; }

	std::map<size_t, std::shared_ptr<Task>>& getLoadedTasks() { return loadedTasks; }

	Loader(std::shared_ptr<AbstractFilter> filter, InputFile::Format inputFormat, int _num_threads, bool storePositions = false);
	
	~Loader() {
		readerQueue.MarkCompleted();
		prefetcherQueue.MarkCompleted();
		intermediateQueue.MarkCompleted();

		for (auto& t : readers) {
			t.join();
		}

		prefetcher.join();
	}

	void configure(const std::string& multipleKmcSamples);
	void initPrefetch();
	void waitForPrefetch() { 
		LOG_VERBOSE << "Waiting for prefetch" << std::endl;
		prefetcherSemaphore.waitForZero();
		LOG_VERBOSE << "Prefetch finished" << std::endl;
	}
	
	void initLoad();
	void waitForLoad() { 
		LOG_VERBOSE << "Waiting for load" << std::endl;
		readerSemaphore.waitForZero(); 
		LOG_VERBOSE << "Load finished" << std::endl;
	
	}

	size_t getBytes() {
		size_t mem = 0;
		for (const auto& col : kmersCollections) {
			mem += col.capacity() * sizeof(kmer_t);
		}
		
		return mem;
	}
	
private:
	
	size_t currentFileId{ 0 };

	InputFile::Format inputFormat;

	int numThreads;

	uint32_t kmerLength;

	bool storePositions;

	std::vector<std::string> kmcFileList;

	std::vector<std::vector<kmer_t>> kmersCollections;

	std::vector<std::vector<uint32_t>> positionsCollections;

	std::thread prefetcher;

	std::vector<std::thread> readers;

	Semaphore prefetcherSemaphore;

	Semaphore readerSemaphore;

	std::mutex outputMutex;

	RegisteringQueue<std::shared_ptr<Task>> readerQueue{ 1 };
	
	RegisteringQueue<std::shared_ptr<Task>> prefetcherQueue{ 1 };

	RegisteringQueue<std::shared_ptr<Task>> intermediateQueue{ 1 };

	std::map<size_t, std::shared_ptr<Task>> loadedTasks;
};