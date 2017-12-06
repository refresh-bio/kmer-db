#include "loader.h"

#include <string>
#include <sstream>
#include <memory>

using namespace std;

/****************************************************************************************************************************************************/
Loader::Loader(std::shared_ptr<IKmerFilter> filter, bool tryMinHash) :
	prefetcherQueue(1), 
	intermediateQueue(1),
	readerQueue(1), 
	currentFileId(0), 
	numThreads(std::thread::hardware_concurrency()),
	filter (filter),
	tryMinHash(tryMinHash) {

	kmersCollections.resize(numThreads);

	// generate preloader thread
	prefetcher = std::thread([this]() {
		while (!this->prefetcherQueue.IsCompleted()) {
			std::shared_ptr<Task> task;
			if (this->prefetcherQueue.Pop(task)) {
				ostringstream oss;
				oss << task->sampleName << " (" << task->fileId + 1 << "/" << kmcFileList.size() << ")...";
				
				task->file = std::make_shared<KmcFileWrapper>(this->filter);
				if (task->file->open(kmcFileList[task->fileId], this->tryMinHash)) {
					intermediateQueue.Push(task);
					oss << "OK";
				}
				else {
					oss << "failed";
				}
				cout << oss.str() << endl;
				prefetcherSemaphore.dec();
			}
		}
	});

	readers.resize(numThreads);
	for (auto& t : readers) {
		t = std::thread([this]() {
			while (!this->readerQueue.IsCompleted()) {
				std::shared_ptr<Task> task;
				if (this->readerQueue.Pop(task)) {
					if (task->file->load(kmersCollections[task->threadId])) {
						task->kmers = &kmersCollections[task->threadId];
						std::unique_lock<std::mutex> lck(outputMutex, std::defer_lock);
						lck.lock();
						loadedTasks[task->fileId] = task;
						lck.unlock();
					}
					readerSemaphore.dec();
				}
			}
		});
	}
}


/****************************************************************************************************************************************************/
void Loader::configure(const std::string& multipleKmcSamples) {
	std::ifstream ifs(multipleKmcSamples);

	string fname;
	while (ifs >> fname) {
		kmcFileList.push_back(fname);
	}

//	sort(kmcFileList.begin(), kmcFileList.end());
//	kmcFileList.erase(unique(kmcFileList.begin(), kmcFileList.end()), kmcFileList.end());

//	kmcFileList.resize(3);
}

/****************************************************************************************************************************************************/
void Loader::initPrefetch() {
	
	for (int tid = 0; tid < numThreads; ++tid) {
		size_t file_id = currentFileId + tid;
		if (file_id < kmcFileList.size()) {

			std::shared_ptr<Task> task = std::make_shared<Task>(file_id, tid, kmcFileList[file_id]);
			
			prefetcherSemaphore.inc();
			prefetcherQueue.Push(task);
		}
	}

	currentFileId = std::min(kmcFileList.size(), currentFileId + numThreads);
}


/****************************************************************************************************************************************************/
void Loader::initLoad() {
	// rewrite stuff from intermediate to reader queue
	std::shared_ptr<Task> task;

	while (!intermediateQueue.IsEmpty()) {
		intermediateQueue.Pop(task);
		readerSemaphore.inc();
		readerQueue.Push(task);
	}
}
