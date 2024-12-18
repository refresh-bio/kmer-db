#include "console.h"
#include "prefix_kmer_db.h"
#include "similarity_calculator.h"
#include "loader_ex.h"
#include "kmer_extract.h"

#include <chrono>
#include <cstdint>
#include <thread>
#include <atomic>

void New2AllConsole::run(const Params& params)
{
	if (params.files.size() != 3) {
		throw usage_error(params.mode);
	}

	LOG_NORMAL("Set of new samples  (from " << InputFile::format2string(params.inputFormat) << ") versus entire database comparison" << endl);

	const std::string& dbFilename = params.files[0];
	const std::string& multipleSamples = params.files[1];
	const std::string & similarityFile = params.files[2];
	
	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(params.numThreads);
	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	std::chrono::duration<double> loadingTime{ 0 }, processingTime{ 0 }, dt{ 0 };

	LOG_NORMAL("Loading k-mer database " << dbFilename << "..." << endl);
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename);
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL("OK (" << dt.count() << " seconds)" << endl << db.printStats() << endl);

	LOG_DEBUG("Creating Loader object..." << endl);
	std::shared_ptr<MinHashFilter> filter(FilterFactory::create(db.getFraction(), db.getStartFraction(), db.getKmerLength()));
	std::shared_ptr<Alphabet> alphabet(AlphabetFactory::instance().create(db.getAlphabetType()));

	LoaderEx loader(filter, alphabet, params.inputFormat, params.numReaderThreads, params.numThreads, params.multisampleFasta);
	loader.configure(multipleSamples);
	LOG_NORMAL(endl);

	std::vector<uint32_t> sims;

	LOG_NORMAL("Processing queries..." << endl);
	auto totalStart = std::chrono::high_resolution_clock::now();

	// create set of buffers for storing similarities
	std::vector<std::vector<uint32_t>> buffers(loader.getOutputBuffersCount());
	RegisteringQueue<int> freeBuffersQueue(1);
	for (size_t i = 0; i < buffers.size(); ++i) {
		freeBuffersQueue.Push(i);
	}

	SynchronizedPriorityQueue<std::shared_ptr<SampleTask>> similarityQueue(params.numThreads);
	std::vector<thread> workers(params.numThreads);
	std::atomic<int> sample_id{ 0 };

	for (int tid = 0; tid < params.numThreads; ++tid) {
		workers[tid] = thread([&db, &loader, &freeBuffersQueue, &similarityQueue, &buffers, &calculator, &sample_id, tid]() {
			int task_id = sample_id.fetch_add(1);
			while (!loader.isCompleted()) {
				std::shared_ptr<SampleTask> task;
				if ((task = loader.popTask(task_id)) && freeBuffersQueue.Pop(task->bufferId2)) {
					LOG_DEBUG("loader queue " << task_id + 1 << " -> (" << task->id + 1 << ", " << task->sampleName << ")" << endl);
					buffers[task->bufferId2].clear();

					// only unique k-mers are needed
					KmerHelper::unique(task->kmers, task->kmersCount);

					calculator.one2all<false>(db, task->kmers, task->kmersCount, buffers[task->bufferId2]);
					similarityQueue.Push(task_id, task);

					LOG_DEBUG("(" << task->id + 1 << ", " << task->sampleName << ") -> similarity queue, tid:" << tid << ", buf:" << task->bufferId2 << endl);
					task_id = sample_id.fetch_add(1);

				}
			}

			similarityQueue.MarkCompleted();
			LOG_DEBUG("similarity thread completed: " << tid << endl);
			});
	}


	// Opening file
	std::ofstream ofs(similarityFile);
	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl;

	// allocate row buffer (10000 for sample name + 100 for each row)
	char* row = new char[10000 + db.getSamplesCount() * 100];
	char* ptr = row;

	ptr += sprintf(ptr, "query-samples,total-kmers,");
	ptr += num2str(db.getSampleKmersCount().data(), db.getSampleKmersCount().size(), ',', ptr);
	*ptr++ = '\n';
	ofs.write(row, ptr - row);

	// Gather results in one thread
	for (int task_id = 0; !similarityQueue.IsCompleted(); ++task_id) {

		std::shared_ptr<SampleTask> task;
		if (similarityQueue.Pop(task_id, task)) {

			if ((task_id + 1) % 10 == 0) {
				LOG_NORMAL("\r" << task_id + 1 << "...                      ");
			}

			LOG_DEBUG("similarity queue -> (" << task_id + 1 << ", " << task->sampleName << "), buf:" << task->bufferId2 << endl);
			auto& buf = buffers[task->bufferId2];

			ptr = row;
			ptr += sprintf(ptr, "%s,%lu,", task->sampleName.c_str(), task->kmersCount);

			if (params.sparseOut) {
				
				std::vector<num_kmers_t> queryKmersCounts(1, task->kmersCount);
				CombinedFilter<num_kmers_t> filter(
					params.metricFilters,
					params.kmerFilter,
					queryKmersCounts,
					db.getSampleKmersCount(),
					db.getKmerLength());
				
				// filter row
				for (int j = 0; j < buf.size(); ++j) {
					if (!filter(buf[j], 0, j)) {
						buf[j] = 0;
					}
				}
				
				ptr += num2str_sparse(buf.data(), buf.size(), ',', ptr);
			}
			else {
				ptr += num2str(buf.data(), buf.size(), ',', ptr);
			}

			freeBuffersQueue.Push(task->bufferId2);
			loader.releaseTask(*task);

			*ptr++ = '\n';
			ofs.write(row, ptr - row);
		}
	}

	delete[] row;

	// make sure all threads have finished
	for (auto& w : workers) {
		w.join();
	}

	auto totalTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);

	LOG_NORMAL(endl << endl << "EXECUTION TIMES" << endl
		<< "Total: " << totalTime.count() << endl);
}
