/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include "loader.h"

#include <string>
#include <sstream>
#include <memory>

using namespace std;

// *****************************************************************************************
//
Loader::Loader(std::shared_ptr<IKmerFilter> filter, InputFile::Format inputFormat, int _num_threads) :
	prefetcherQueue(1), 
	intermediateQueue(1),
	readerQueue(1), 
	currentFileId(0), 
	numThreads(_num_threads > 0 ? _num_threads : std::thread::hardware_concurrency()),
	inputFormat(inputFormat)
	{
	kmersCollections.resize(numThreads);
	
	// generate preloader thread
	prefetcher = std::thread([this, filter]() {
		while (!this->prefetcherQueue.IsCompleted()) {
			std::shared_ptr<Task> task;
			if (this->prefetcherQueue.Pop(task)) {
//				cout << "\r" << std::string(task->sampleName.size() + 100, ' ') << "\r";
				ostringstream oss;
				cout << "\r" << task->sampleName << " (" << task->fileId + 1 << "/" << kmcFileList.size() << ")...                      ";
				fflush(stdout);

				if (this->inputFormat == InputFile::KMC) {
					task->file = std::make_shared<KmcInputFile>(filter->clone());
				}
				else if (this->inputFormat == InputFile::MINHASH) {
					task->file = std::make_shared<MihashedInputFile>(filter->clone());
				}
				else {
					task->file = std::make_shared<GenomeInputFile>(filter->clone());
				}

				if (task->file->open(kmcFileList[task->fileId])) {
					intermediateQueue.Push(task);
				}
				else {
					cout << "failed" << endl;
				}
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
					if (task->file->load(kmersCollections[task->threadId], task->kmerLength, task->fraction)) {
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

// *****************************************************************************************
//
void Loader::configure(const std::string& multipleKmcSamples) {
	std::ifstream ifs(multipleKmcSamples);

	string fname;
	while (ifs >> fname) {
		kmcFileList.push_back(fname);
	}
}

// *****************************************************************************************
//
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

// *****************************************************************************************
//
void Loader::initLoad() {
	// rewrite stuff from intermediate to reader queue
	std::shared_ptr<Task> task;

	while (!intermediateQueue.IsEmpty()) {
		intermediateQueue.Pop(task);
		readerSemaphore.inc();
		readerQueue.Push(task);
	}
}
