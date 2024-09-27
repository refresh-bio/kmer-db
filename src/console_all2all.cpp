#include "console.h"
#include "similarity_calculator.h"

#include <chrono>
#include <cstdint>

void All2AllConsole::run(const Params& params) {

	if (params.files.size() != 2) {
		throw usage_error(params.mode);
	}
	
	LOG_NORMAL << "All versus all comparison" << endl;

	const std::string& dbFilename = params.files[0];
	const std::string& similarityFile = params.files[1];

	std::ifstream dbFile(dbFilename, std::ios::binary);
	std::ofstream ofs(similarityFile);
	PrefixKmerDb* db = new PrefixKmerDb(params.numThreads);
	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	std::chrono::duration<double> dt{ 0 };
	LOG_NORMAL << "Loading k-mer database " << dbFilename << "..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile || !db->deserialize(dbFile, AbstractKmerDb::DeserializationMode::SkipHashtables)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename);
	}
	dt = std::chrono::high_resolution_clock::now() - start;

	LOG_NORMAL << "Calculating matrix of common k-mers..." << endl;
	start = std::chrono::high_resolution_clock::now();
	LowerTriangularMatrix<uint32_t> matrix;
	calculator.all2all(*db, matrix);
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

	if (params.sparseOut) {

		CombinedFilter<uint32_t> filter(
			params.metricFilters,
			params.kmerFilter,
			db->getSampleKmersCount(),
			db->getSampleKmersCount(),
			params.kmerLength);

		matrix.compact(filter);
	}

	for (size_t sid = 0; sid < db->getSamplesCount(); ++sid) {
		ptr = row;
		ptr += sprintf(ptr, "%s,%lu,", db->getSampleNames()[sid].c_str(), (unsigned long)db->getSampleKmersCount()[sid]);

		if (params.sparseOut) {
			ptr += matrix.saveRowSparse(sid, ptr);
		}
		else {
			ptr += matrix.saveRow(sid, ptr);
		}

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