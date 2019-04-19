#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "input_file.h"
#include "queue.h"
#include "filter.h"
#include "loader_tasks.h"

#include <vector>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>



// *****************************************************************************************
//
class LoaderEx {
public:

	LoaderEx(
		std::shared_ptr<AbstractFilter> filter, 
		InputFile::Format inputFormat, 
		int numThreads,
		bool multisampleFasta,
		bool storePositions = false);

	~LoaderEx();

	int configure(const std::string& multipleKmcSamples);

	std::shared_ptr<SampleTask> popTask(int fileId) {
		std::shared_ptr<SampleTask> task;
		queues.output.Pop(fileId, task);
		LOG_DEBUG << "output queue -> (" << fileId + 1 << ")" << std::endl << std::flush;
		return task;
	}

	void releaseTask(SampleTask& t) {
		if (--bufferRefCounters[t.bufferId] == 0) {
			queues.freeBuffers.Push(t.bufferId);
			LOG_DEBUG << "Released readers buffer: " << t.fileId + 1 << std::endl << std::flush;
		}
	}
	
	size_t getBytes() {
		size_t mem = 0;
		for (const auto& col : kmersCollections) {
			mem += col.capacity() * sizeof(kmer_t);
		}

		return mem;
	}

	bool isCompleted() {
		return queues.output.IsCompleted();
	}

private:

	InputFile::Format inputFormat;

	int numThreads;

	bool multisampleFasta;

	bool storePositions;

	uint32_t kmerLength;

	std::thread prefetcher;

	std::vector<std::string> fileNames;

	std::vector<std::thread> readers;

	std::vector<std::vector<kmer_t>> kmersCollections;

	std::vector<std::vector<uint32_t>> positionsCollections;

	std::vector<int> bufferRefCounters;

	std::atomic<sample_id_t> multisampleCounter{ 0 };

	struct {
		RegisteringQueue<std::shared_ptr<InputTask>> input{ 1 };

		RegisteringQueue<std::shared_ptr<InputTask>> readers{ 1 };

		RegisteringQueue<int> freeBuffers{ 1 };

		SynchronizedPriorityQueue<std::shared_ptr<SampleTask>> output{ 1 };
	} queues;

};