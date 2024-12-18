#include "console.h"
#include "minhashed_input_file.h"
#include "loader_ex.h"
#include "kmer_extract.h"

void MinhashConsole::run(const Params& params) {
	
	if (params.files.size() != 1) {
		throw usage_error(params.mode);
	}

	LOG_NORMAL("Minhashing samples..." << endl);

	const std::string& multipleKmcSamples = params.files[0];
	std::chrono::duration<double> loadingTime{ 0 }, processingTime{ 0 };

	LOG_DEBUG("Creating Loader object..." << endl);

	std::shared_ptr<MinHashFilter> filter(FilterFactory::create(params.fraction, 0, params.kmerLength));

	LoaderEx loader(filter, params.alphabet, params.inputFormat, params.numReaderThreads, params.numThreads, params.multisampleFasta);
	loader.configure(multipleKmcSamples);

	LOG_DEBUG("Starting loop..." << endl);
	auto totalStart = std::chrono::high_resolution_clock::now();
	for (int i = 0; !loader.isCompleted(); ++i) {
		auto partialTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart);
		LOG_VERBOSE("Processing time: " << partialTime.count() << ", loader buffers: " << (loader.getBytes() >> 20) << " MB" << endl);

		auto task = loader.popTask(i);

		if (task) {
			auto start = std::chrono::high_resolution_clock::now();

			MihashedInputFile file;

			// postprocess k-mers if neccessary
			if (params.inputFormat == InputFile::Format::GENOME) {
				KmerHelper::sortAndUnique(task->kmers, task->kmersCount, params.numThreads);
			}
			else if (params.inputFormat == InputFile::Format::KMC) {
				KmerHelper::sort(task->kmers, task->kmersCount, params.numThreads);
			}

			file.store(task->filePath, task->kmers, task->kmersCount, task->kmerLength, params.fraction);

			processingTime += std::chrono::high_resolution_clock::now() - start;
			loader.releaseTask(*task);
		}
	}
}