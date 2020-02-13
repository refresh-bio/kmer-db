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
	multisampleFasta(multisampleFasta),
	storePositions(storePositions)
{
	readers.resize(numThreads);
	int buffersCount = std::max(numThreads, 8);

	// configure queues
	queues.input.Restart(1);
	queues.readers.Restart(1, buffersCount);
	queues.output.Restart(numThreads);
	queues.freeBuffers.Restart(1);
	
	kmersCollections.resize(buffersCount);
	positionsCollections.resize(buffersCount);
	bufferRefCounters.resize(buffersCount);
	
	for (size_t id = 0; id < kmersCollections.size(); ++id) {
		queues.freeBuffers.Push((int)id);
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

		queues.readers.MarkCompleted();
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

					LOG_DEBUG << "readers queue -> (" << sampleTask->fileId + 1 << "), tid: " << tid << ", buf: " << bufferId << endl << std::flush;
					if ((sampleTask->fileId + 1) % 10 == 0) {
						cout << "\r" << sampleTask->fileId + 1 << "/" << fileNames.size() << "...                      " << std::flush;
					}

					bool ok = false;

					if (this->multisampleFasta) {
						auto genomicFile = std::dynamic_pointer_cast<GenomeInputFile>(inputTask->file);
						int count = genomicFile->loadMultiple(kmersCollections[bufferId], positionsCollections[bufferId], 
							sampleTask, this->multisampleCounter, queues.output);
						bufferRefCounters[bufferId] += count;
						ok = (count > 0);
					}
					else {
						ok = inputTask->file->load(
							kmersCollections[bufferId], positionsCollections[bufferId],
							sampleTask->kmers, sampleTask->kmersCount, sampleTask->kmerLength, sampleTask->fraction);
						
						if (ok) {
							++bufferRefCounters[bufferId];
							queues.output.Push((int)sampleTask->fileId, sampleTask);
						}
					}

					if (ok) {
						LOG_VERBOSE << "File loaded successfully: " << sampleTask->fileId + 1 << endl << std::flush;
					} else {
						cout << "File load failed: " << sampleTask->fileId + 1 << endl << std::flush;
					}
				}
			}

			queues.output.MarkCompleted();
		});
	}
}

// *****************************************************************************************
//
LoaderEx::~LoaderEx() {
	queues.input.MarkCompleted();
	queues.readers.MarkCompleted();
	queues.output.MarkCompleted();	
	queues.freeBuffers.MarkCompleted();

	for (auto& t : readers) {
		t.join();
	}

	prefetcher.join();
}

// *****************************************************************************************
//
int LoaderEx::configure(const std::string& multipleSamples) {
	std::ifstream ifs(multipleSamples);

	string fname;
	
	while (ifs >> fname) {
		fileNames.push_back(fname);
	}

	for (size_t i = 0; i < fileNames.size(); ++i) {
		queues.input.Push(std::make_shared<InputTask>(i, fileNames[i]));
	}
	
	queues.input.MarkCompleted();
	return (int)fileNames.size();
}


