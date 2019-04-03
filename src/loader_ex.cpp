/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "loader_ex.h"

#include <string>
#include <sstream>
#include <memory>

using namespace std;

// *****************************************************************************************
//
LoaderEx::LoaderEx(
	std::shared_ptr<AbstractFilter> filter, 
	InputFile::Format inputFormat, 
	int _num_threads, 
	bool storePositions) :

	inputFormat(inputFormat),
	numThreads(_num_threads),
	storePositions(storePositions)

{
	readers.resize(numThreads);
	int buffersCount = std::max(numThreads, 8);

	queues.readers.Restart(1, buffersCount);
	kmersCollections.resize(buffersCount);
	positionsCollections.resize(buffersCount);
	bufferRefCounters.resize(buffersCount);
	
	for (int id = 0; id < kmersCollections.size(); ++id) {
		queues.freeBuffers.Push(id);
		kmersCollections[id].reserve(10000000);
		if (storePositions) {
			positionsCollections[id].reserve(10000000);
		}
	}


	// generate preloader thread
	prefetcher = std::thread([this, filter]() {
		while (!this->queues.input.IsCompleted()) {
			std::shared_ptr<TaskEx> task;
			if (this->queues.input.Pop(task)) {
				LOG_DEBUG << "input queue -> (" << task->fileId + 1 << ")" << endl << std::flush;
				
				if (this->inputFormat == InputFile::KMC) {
					task->file = std::make_shared<KmcInputFile>(filter->clone());
				}
				else if (this->inputFormat == InputFile::MINHASH) {
					task->file = std::make_shared<MihashedInputFile>(filter->clone());
				}
				else {
					task->file = std::make_shared<GenomeInputFile>(filter->clone(), this->storePositions);
				}

				if (task->file->open(task->filePath)) {
					queues.readers.Push(task);
					LOG_DEBUG << "(" << task->fileId + 1  << ") -> readers queue " << endl << std::flush;
				}
				else {
					cout << "failed:" << task->sampleName << endl << std::flush;
				}
			}
		}
	});

	for (int tid = 0; tid < numThreads; ++tid) {
		readers[tid] = std::thread([this, tid]() {
			while (!this->queues.readers.IsCompleted()) {
				std::shared_ptr<TaskEx> task;
				int bufferId;

				if (this->queues.freeBuffers.Pop(bufferId) && this->queues.readers.Pop(task)) {
					
					task->bufferId = bufferId;
					
					LOG_DEBUG << "readers queue -> (" << task->fileId + 1 << "), tid: " << tid << ", buf: " << bufferId  << endl << std::flush;
					if ((task->fileId + 1) % 10 == 0) {
						cout << "\r" << task->fileId + 1 << "/" << numFiles << "...                      " << std::flush;
					}
					
					if (task->file->load(kmersCollections[bufferId], positionsCollections[bufferId], task->kmerLength, task->fraction)) {
						LOG_VERBOSE << "Sample loaded successfully: " << task->fileId + 1 << endl << std::flush;
						task->kmers = kmersCollections[bufferId].data();
						task->kmersCount = kmersCollections[bufferId].size();
						++bufferRefCounters[bufferId];
						queues.output.Push(task->fileId, task);

					}
					else {
						cout << "Sample load failed: " << task->fileId + 1 << endl << std::flush;
					}
				}
			}
		});
	}
}

// *****************************************************************************************
//
LoaderEx::~LoaderEx() {
	queues.input.MarkCompleted();
	queues.readers.MarkCompleted();
	queues.output.MarkCompleted();

	for (auto& t : readers) {
		t.join();
	}

	prefetcher.join();
}

// *****************************************************************************************
//
int LoaderEx::configure(const std::string& multipleKmcSamples) {
	std::ifstream ifs(multipleKmcSamples);

	string fname;
	int fid = 0;
	while (ifs >> fname) {
		queues.input.Push(std::make_shared<TaskEx>(fid, fname));
		++fid;
	}

	this->numFiles = fid;

	return fid;
}


