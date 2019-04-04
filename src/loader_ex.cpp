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
	int numThreads, 
	bool multisampleFasta,
	bool storePositions) :

	inputFormat(inputFormat),
	numThreads(numThreads),
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
			std::shared_ptr<InputTask> task;
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
					cout << "failed:" << task->filePath << endl << std::flush;
				}
			}
		}
	});

	for (int tid = 0; tid < numThreads; ++tid) {
		readers[tid] = std::thread([this, tid]() {
			while (!this->queues.readers.IsCompleted()) {
				std::shared_ptr<InputTask> inputTask;
				int bufferId;

				if (this->queues.freeBuffers.Pop(bufferId) && this->queues.readers.Pop(inputTask)) {
					
					auto sampleTask = make_shared<SampleTask>(
						inputTask->fileId,
						inputTask->filePath,
						InputFile::removePathFromFile(inputTask->filePath),
						bufferId);
					
					LOG_DEBUG << "readers queue -> (" << sampleTask->fileId + 1 << "), tid: " << tid << ", buf: " << bufferId  << endl << std::flush;
					if ((sampleTask->fileId + 1) % 10 == 0) {
						cout << "\r" << sampleTask->fileId + 1 << "/" << numFiles << "...                      " << std::flush;
					}
				
					if (inputTask->file->load(
						kmersCollections[bufferId], positionsCollections[bufferId], 
						sampleTask->kmers, sampleTask->kmersCount, sampleTask->kmerLength, sampleTask->fraction)) {
						
						LOG_VERBOSE << "Sample loaded successfully: " << sampleTask->fileId + 1 << endl << std::flush;
						++bufferRefCounters[bufferId];
						queues.output.Push(sampleTask->fileId, sampleTask);
					}
					else {
						cout << "Sample load failed: " << sampleTask->fileId + 1 << endl << std::flush;
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
		queues.input.Push(std::make_shared<InputTask>(fid, fname));
		++fid;
	}

	this->numFiles = fid;

	return fid;
}


