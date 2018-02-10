#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include "kmc_file_wrapper.h"
#include "queue.h"
#include "filter.h"

#include <vector>
#include <memory>
#include <map>
#include <fstream>

// *****************************************************************************************
//
struct Task {
	size_t fileId;
	size_t threadId;
	std::string filePath;
	std::string sampleName;
	std::shared_ptr<KmcFileWrapper> file;
	std::vector<kmer_t>* kmers;
	uint32_t kmerLength;
	double fraction;

	// *****************************************************************************************
	//
	Task(size_t fileId, size_t threadId, const std::string& filePath) :
		fileId(fileId), threadId(threadId), filePath(filePath), file(nullptr), kmers(nullptr) {
	
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
	
	Loader(std::shared_ptr<IKmerFilter> filter, bool tryMinHash, int _num_threads);
	// *****************************************************************************************
	//
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
	void waitForPrefetch() { prefetcherSemaphore.waitForZero();  }
	void initLoad();
	void waitForLoad() { readerSemaphore.waitForZero();  }

	std::map<size_t, std::shared_ptr<Task>>& getLoadedTasks() { return loadedTasks;  }

private:
	
	bool useMinhash;

	size_t currentFileId;

	int numThreads;

	std::vector<std::string> kmcFileList;

	std::vector<std::vector<kmer_t>> kmersCollections;

	std::thread prefetcher;

	std::vector<std::thread> readers;

	Semaphore prefetcherSemaphore;

	Semaphore readerSemaphore;

	std::mutex outputMutex;

	CRegisteringQueue<std::shared_ptr<Task>> readerQueue;
	
	CRegisteringQueue<std::shared_ptr<Task>> prefetcherQueue;

	CRegisteringQueue<std::shared_ptr<Task>> intermediateQueue;

	std::map<size_t, std::shared_ptr<Task>> loadedTasks;
};