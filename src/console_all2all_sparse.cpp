#include "console.h"
#include "prefix_kmer_db.h"
#include "similarity_calculator.h"

#include <chrono>
#include <cstdint>
#include <algorithm>

// *****************************************************************************************
//
void All2AllSparseConsole::run(const Params& params) {

	
	if (params.files.size() != 2) {
		throw usage_error(params.mode);
	}

	LOG_NORMAL << "All versus all comparison (sparse computation)" << endl;

	const std::string& dbFilename = params.files[0];
	const std::string& similarityFile = params.files[1];
	
	//uint32_t below = (uint32_t)lrint(params.below);
	//uint32_t above = (uint32_t)std::max(0l, lrint(params.above));

	std::ifstream dbFile(dbFilename, std::ios::binary);
	std::ofstream ofs(similarityFile, std::ios::binary);
	PrefixKmerDb* db = new PrefixKmerDb(params.numThreads);
	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	std::chrono::duration<double> dt{ 0 };
	LOG_NORMAL << "Loading k-mer database " << dbFilename << "..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db->deserialize(dbFile, AbstractKmerDb::DeserializationMode::SkipHashtables)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename);
	}
	dt = std::chrono::high_resolution_clock::now() - start;

	LOG_NORMAL << "Calculating matrix of common k-mers...";
	start = std::chrono::high_resolution_clock::now();
	SparseMatrix<uint32_t> matrix;
	calculator.all2all_sp(*db, matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;

	LOG_NORMAL << "Storing matrix of common k-mers in " << similarityFile << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db->getKmerLength() << " fraction: " << db->getFraction() << " ,db-samples ,";
	std::copy(db->getSampleNames().cbegin(), db->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl;

	// allocate row buffer (10000 for sample name + 100 for each row)
	char* row = new char[10000 + db->getSamplesCount() * 100];
	char* ptr = row;

	ptr += sprintf(ptr, "query-samples,total-kmers,");
	ptr += num2str(db->getSampleKmersCount().data(), db->getSampleKmersCount().size(), ',', ptr);
	*ptr++ = '\n';
	ofs.write(row, ptr - row);

	for (size_t sid = 0; sid < db->getSamplesCount(); ++sid) {
		ptr = row;
		ptr += sprintf(ptr, "%s,%lu,", db->getSampleNames()[sid].c_str(), db->getSampleKmersCount()[sid]);

		ptr += matrix.saveRowSparse(sid, ptr);

		*ptr++ = '\n';
		ofs.write(row, ptr - row);
	}

	delete[] row;

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;

	LOG_NORMAL << "Releasing memory...";
	start = std::chrono::high_resolution_clock::now();
	delete db;
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;
}