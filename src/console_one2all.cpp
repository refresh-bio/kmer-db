#include "console.h"
#include "prefix_kmer_db.h"
#include "similarity_calculator.h"
#include "loader_ex.h"
#include "kmer_extract.h"

#include "input_file_factory.h"

#include <chrono>
#include <cstdint>

void One2AllConsole::run(const Params& params) {

	if (params.files.size() != 3) {
		throw usage_error(params.mode);
	}
	
	LOG_NORMAL("One new sample  (from " << InputFile::format2string(params.inputFormat) << ") versus entire database comparison" << endl);
	
	const std::string& dbFilename = params.files[0];
	const std::string& sampleFasta = params.files[1];
	const std::string& similarityFile = params.files[2];
	
	//uint32_t below = (uint32_t)lrint(params.below);
	//uint32_t above = (uint32_t)std::max(0l, lrint(params.above));

	std::ifstream dbFile(dbFilename, std::ios::binary);
	PrefixKmerDb db(params.numThreads);
	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	std::chrono::duration<double> dt{ 0 };

	LOG_NORMAL("Loading k-mer database " << dbFilename << ":" << endl);
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db.deserialize(dbFile)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename);
	}
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL("OK (" << dt.count() << " seconds)" << endl << db.printStats() << endl);

	LOG_NORMAL("Loading sample kmers...");

	start = std::chrono::high_resolution_clock::now();

	std::vector<kmer_t> kmersBuffer;
	std::vector<uint32_t> positions;
	uint32_t kmerLength;
	
	std::shared_ptr<MinHashFilter> filter(FilterFactory::create(db.getFraction(), db.getStartFraction(), db.getKmerLength()));
	std::shared_ptr<Alphabet> alphabet(AlphabetFactory::instance().create(db.getAlphabetType()));

	std::shared_ptr<InputFile> file(InputFileFactory::create(params.inputFormat, filter, alphabet));

	double dummy;
	size_t queryKmersCount;
	kmer_t* queryKmers;

	if (!file->open(sampleFasta) || !file->load(kmersBuffer, positions, queryKmers, queryKmersCount, kmerLength, dummy)) {
		throw runtime_error("Cannot open sample file: " + sampleFasta);
	}

	// postprocess k-mers if neccessary
	if (params.inputFormat == InputFile::Format::GENOME) {
		KmerHelper::sortAndUnique(queryKmers, queryKmersCount, params.numThreads);
	}

	if (kmerLength != db.getKmerLength()) {
		throw runtime_error("Sample and database k-mer length differ");
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL("OK (" << dt.count() << " seconds)" << endl
		<< "Number of k-mers: " << queryKmersCount << endl
		<< "Minhash fraction: " << db.getFraction() << endl);

	LOG_NORMAL("Calculating similarity vector...");
	start = std::chrono::high_resolution_clock::now();
	std::vector<uint32_t> sims;
	calculator.one2all(db, queryKmers, queryKmersCount, sims);
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL("OK (" << dt.count() << " seconds)" << endl);

	LOG_NORMAL("Storing similarity vector in " << similarityFile << "...");
	std::ofstream ofs(similarityFile);

	ofs << "kmer-length: " << db.getKmerLength() << " fraction: " << db.getFraction() << " ,db-samples ,";
	std::copy(db.getSampleNames().cbegin(), db.getSampleNames().cend(), ostream_iterator<string>(ofs, ","));

	ofs << endl << "query-samples,total-kmers,";
	std::copy(db.getSampleKmersCount().cbegin(), db.getSampleKmersCount().cend(), ostream_iterator<size_t>(ofs, ","));
	ofs << endl << sampleFasta << "," << queryKmersCount << ",";
	std::copy(sims.begin(), sims.end(), ostream_iterator<uint32_t>(ofs, ","));

	ofs.close();
	LOG_NORMAL("OK" << endl);
}

