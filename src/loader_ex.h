#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "kmc_file_wrapper.h"
#include "queue.h"
#include "filter.h"
#include "loader.h"

#include <vector>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>



// *****************************************************************************************
//
struct TaskEx {
	size_t fileId;
	std::string filePath;
	std::string sampleName;
	std::shared_ptr<InputFile> file;
	kmer_t *kmers;
	size_t kmersCount;
	uint32_t kmerLength;
	double fraction;
	int bufferId;

	// *****************************************************************************************
	//
	TaskEx(size_t fileId, const std::string& filePath) :
		fileId(fileId), filePath(filePath), file(nullptr), kmers(nullptr), kmersCount(0) {

		size_t pos = filePath.find_last_of("/\\");
		if (pos != string::npos) {
			sampleName = filePath.substr(pos + 1);
		}
		else {
			sampleName = filePath;
		}
	}

	bool operator<(const TaskEx& other) {
		return this->fileId > other.fileId;
	}
};

// *****************************************************************************************
//
class LoaderEx {
public:

	LoaderEx(
		std::shared_ptr<AbstractFilter> filter, 
		InputFile::Format inputFormat, 
		int _num_threads, 
		bool storePositions = false);

	~LoaderEx();

	int configure(const std::string& multipleKmcSamples);

	std::shared_ptr<TaskEx> popTask(int fileId) {
		std::shared_ptr<TaskEx> task;
		queues.output.Pop(fileId, task);
		LOG_DEBUG << "output queue -> (" << fileId + 1 << ")" << std::endl << std::flush;
		return task;
	}

	void releaseTask(TaskEx& t) {
		queues.freeBuffers.Push(t.bufferId);
		LOG_DEBUG << "Released readers buffer: " << t.fileId + 1 << std::endl << std::flush;
	}
	
	size_t getBytes() {
		size_t mem = 0;
		for (const auto& col : kmersCollections) {
			mem += col.capacity() * sizeof(kmer_t);
		}

		return mem;
	}

private:

	InputFile::Format inputFormat;

	int numThreads;

	bool storePositions;

	int numFiles;

	uint32_t kmerLength;

	std::thread prefetcher;

	std::vector<std::thread> readers;


	std::vector<std::vector<kmer_t>> kmersCollections;

	std::vector<std::vector<uint32_t>> positionsCollections;

	struct {
		RegisteringQueue<std::shared_ptr<TaskEx>> input{ 1 };

		RegisteringQueue<std::shared_ptr<TaskEx>> readers{ 1 };

		RegisteringQueue<int> freeBuffers{ 1 };

		SynchronizedPriorityQueue<std::shared_ptr<TaskEx>> output{ 1 };
	} queues;

};