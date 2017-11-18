#include "loader.h"

#include <string>
#include <sstream>

using namespace std;

/****************************************************************************************************************************************************/
Loader::Loader(std::shared_ptr<IKmerFilter> filter) :
	prefetcherQueue(1), 
	intermediateQueue(1),
	readerQueue(1), 
	currentFileId(0), 
	numThreads(std::thread::hardware_concurrency()),
	filter (filter) {

	kmersCollections.resize(numThreads);

	// generate preloader thread
	prefetcher = std::thread([this]() {
		while (!this->prefetcherQueue.IsCompleted()) {
			std::shared_ptr<Task> task;
			if (this->prefetcherQueue.Pop(task)) {
				task->file = std::make_shared<CKMCFile>();
				ostringstream oss;
				oss << task->sampleName << " (" << task->fileId + 1 << "/" << kmcFileList.size() << ")...";
				if (task->file->OpenForListing(kmcFileList[task->fileId])) {
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
					if (loadKmers(*(task->file), kmersCollections[task->threadId])) {
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

	sort(kmcFileList.begin(), kmcFileList.end());
	kmcFileList.erase(unique(kmcFileList.begin(), kmcFileList.end()), kmcFileList.end());

	//kmcFileList.resize(100);
}

/****************************************************************************************************************************************************/
void Loader::initPrefetch() {
	
	for (int tid = 0; tid < numThreads; ++tid) {
		size_t file_id = currentFileId + tid;
		if (file_id < kmcFileList.size()) {

			std::shared_ptr<Task> task = std::make_shared<Task>(file_id, tid, kmcFileList[file_id]);
			
			size_t pos = task->sampleName.find_last_of("/\\");
			if (pos != string::npos) {
				task->sampleName = task->sampleName.substr(pos + 1);
			}

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

/****************************************************************************************************************************************************/
bool Loader::loadKmers(const string &filename, std::vector<kmer_t>& kmers) {
	CKMCFile kmc_file;
	
	if (!kmc_file.OpenForListing(filename))
		return false;

	return loadKmers(kmc_file, kmers);
}

/****************************************************************************************************************************************************/
bool Loader::loadKmers(CKMCFile& kmc_file, std::vector<kmer_t>& kmers) {
	uint32_t counter;
	
	uint32 _kmer_length;
	uint32 _mode;
	uint32 _counter_size;
	uint32 _lut_prefix_length;
	uint32 _signature_len;
	uint32 _min_count;
	uint64 _max_count;
	uint64 _total_kmers;

	kmc_file.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	CKmerAPI kmer(_kmer_length);

	uint64_t u_kmer;
	vector<uint64> tmp;

	// Wczytuje wszystkie k-mery z pliku do wektora, zeby pozniej moc robic prefetcha
	kmers.clear();

	while (!kmc_file.Eof())
	{
		if (!kmc_file.ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		if ((*filter)(u_kmer)) {
			kmers.push_back(u_kmer);
		}
	}

	return kmc_file.Close();
}
