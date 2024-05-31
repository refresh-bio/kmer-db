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
	int suggestedNumThreads, 
	int numConsumers,
	bool multisampleFasta,
	bool storePositions) :

	inputFormat(inputFormat),
	numThreads(multisampleFasta ? 1 : suggestedNumThreads), // use only one reader thread for multisample fasta files
	multisampleFasta(multisampleFasta),
	storePositions(storePositions)
{
	readers.resize(numThreads);
	
	int outputBuffersCount = numThreads + numConsumers;
	
	// configure queues
	queues.input.Restart(1);
	queues.readers.Restart(1, outputBuffersCount);
	queues.output.Restart(numThreads);
	queues.freeBuffers.Restart(1);
	
	kmersCollections.resize(outputBuffersCount);
	positionsCollections.resize(outputBuffersCount);
	bufferRefCounters.resize(outputBuffersCount);
	
	for (size_t id = 0; id < kmersCollections.size(); ++id) {
		queues.freeBuffers.Push((int)id);
	}

	// run prefetcher thread
	prefetcher = std::thread(&LoaderEx::prefetcherJob, this, filter);

	// run loader threads
	if (multisampleFasta) {
		// multisample fasta - single reader thread
		readers[0] = std::thread(&LoaderEx::multifastaReaderJob, this);
	} else {
		// collection of fasta files - several reader threads
		for (int tid = 0; tid < numThreads; ++tid) {
			readers[tid] = std::thread(&LoaderEx::readerJob, this, tid);
		}
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

	if (!ifs) {
		throw std::runtime_error("Unable to open input file " + multisampleFasta);
	}

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

// *****************************************************************************************
//
void LoaderEx::prefetcherJob(std::shared_ptr<AbstractFilter> filter) {
	while (!this->queues.input.IsCompleted()) {
		std::shared_ptr<InputTask> task;
		
		if (this->queues.input.Pop(task)) {
			LOG_DEBUG << "input queue -> (file " << task->fileId + 1 << ")" << endl ;

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
				LOG_DEBUG << "(file " << task->fileId + 1 << ", " << task->filePath << ") -> readers queue " << endl ;
			}
			else {
				LOG_NORMAL << "failed:" << task->filePath << endl;
			}
		}
	}

	queues.readers.MarkCompleted();
	LOG_DEBUG << "reader thread completed" << endl ;
}

// *****************************************************************************************
//
void LoaderEx::readerJob(int tid) {
	
	while (!this->queues.readers.IsCompleted()) {
		std::shared_ptr<InputTask> inputTask;
		int bufferId = 0;
		bool ok = false;

		// get buffer and input task
		if (this->queues.freeBuffers.Pop(bufferId) && this->queues.readers.Pop(inputTask)) {

			LOG_DEBUG << "readers queue -> (file " << inputTask->fileId + 1 << "), tid: " << tid << endl;

			auto sampleTask = make_shared<SampleTask>(
				inputTask->fileId,
				inputTask->filePath,
				InputFile::removePathFromFile(inputTask->filePath),
				bufferId);

			ok = inputTask->file->load(
				kmersCollections[bufferId], positionsCollections[bufferId],
				sampleTask->kmers, sampleTask->kmersCount, sampleTask->kmerLength, sampleTask->fraction, false);

			if (ok) {
				++bufferRefCounters[bufferId];
				queues.output.Push((int)sampleTask->id, sampleTask);
				LOG_DEBUG << "(sample " << sampleTask->id + 1 << ", " << sampleTask->sampleName <<") -> loader output queue, buf: " << bufferId << std::endl ;
				LOG_VERBOSE << "File loaded successfully: " << inputTask->fileId + 1 << endl ;
			}
			else {
				LOG_NORMAL << "File load failed: " << inputTask->fileId + 1 << endl << std::flush;
			}
		}
	}

	queues.output.MarkCompleted();
	LOG_DEBUG << "loader thread completed: " << tid << endl ;
}


// *****************************************************************************************
//
void LoaderEx::multifastaReaderJob() {

	size_t sample_id = 0;

	while (!this->queues.readers.IsCompleted()) {
		std::shared_ptr<InputTask> inputTask;
		int bufferId = 0;
		int count = 0;

		// wait for input task
		if (!this->queues.readers.Pop(inputTask)) {
			continue;
		}

		auto genomicFile = std::dynamic_pointer_cast<GenomeInputFile>(inputTask->file);

		// initialize multifasta file
		if (genomicFile->initMultiFasta()) {

			LOG_DEBUG << "multifasta initialized: " << inputTask->fileId + 1 << endl ;

			while (true) {
				LOG_DEBUG << "wait for buf for sample " << sample_id + 1 << endl;
				this->queues.freeBuffers.Pop(bufferId); // wait for free buffer
				LOG_DEBUG << "acquired buf " << bufferId << " for sample " << sample_id + 1 << endl;
			
				auto sampleTask = make_shared<SampleTask>(
					sample_id,
					inputTask->filePath,
					"",
					bufferId);

				bool ok = genomicFile->loadNext(
					kmersCollections[bufferId], positionsCollections[bufferId],
					sampleTask->kmers, sampleTask->kmersCount, sampleTask->kmerLength, sampleTask->fraction, false, sampleTask->sampleName);

				++sample_id;
				++bufferRefCounters[bufferId];
				queues.output.Push((int)sampleTask->id, sampleTask);
				++count;

				LOG_DEBUG << "(sample " << sampleTask->id + 1 << ") -> output queue, buf: " << bufferId << std::endl;

				// no more samples
				if (!ok) {
					break;
				}
			}
		}

		if (count > 0) {
			LOG_VERBOSE << "File loaded successfully: " << inputTask->fileId + 1 << endl ;
		}
		else {
			LOG_NORMAL << "File load failed: " << inputTask->fileId + 1 << endl << std::flush;
		}
	}
	queues.output.MarkCompleted();

	LOG_DEBUG << "output queue: mark completed" << endl ;
}

