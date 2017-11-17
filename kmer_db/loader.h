#pragma once

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "queue.h"

#include <vector>
#include <memory>
#include <map>


struct Task {
	size_t fileId;
	size_t threadId;
	std::string sampleName;
	std::shared_ptr<CKMCFile> file;
	std::vector<kmer_t>* kmers;

	Task(size_t fileId, size_t threadId, const std::string& sampleName) :
		fileId(fileId), threadId(threadId), sampleName(sampleName), file(nullptr), kmers(nullptr) {}
};


class Loader {
public:
	Loader();
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

	bool loadKmers(CKMCFile& file, std::vector<kmer_t>& kmers);
	bool loadKmers(const string &filename, std::vector<kmer_t>& kmers);

	std::map<size_t, std::shared_ptr<Task>>& getLoadedTasks() { return loadedTasks;  }

private:
	
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