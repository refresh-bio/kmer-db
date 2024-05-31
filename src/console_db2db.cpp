#include "console.h"
#include "prefix_kmer_db.h"
#include "similarity_calculator.h"
#include "loader_ex.h"
#include "kmer_extract.h"

#include <chrono>
#include <cstdint>


void Db2DbConsole::run(const Params& params)
{
	if (params.files.size() != 3) {
		throw usage_error(params.mode);
	}
	
	const std::string& dbFilename1 = params.files[0];
	const std::string& dbFilename2 = params.files[1];
	const std::string& similarityFilename = params.files[2];
	
	std::ifstream dbFile1(dbFilename1, std::ios::binary);
	PrefixKmerDb* db1 = new PrefixKmerDb(params.numThreads);

	std::ifstream dbFile2(dbFilename2, std::ios::binary);
	PrefixKmerDb* db2 = new PrefixKmerDb(params.numThreads);

	std::ofstream ofs(similarityFilename, std::ios::binary);

	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	std::chrono::duration<double> loadingTime{ 0 }, processingTime{ 0 }, dt{ 0 };

	// !!! TODO: zrobiæ równoleg³y odczyt obu baz
	LOG_NORMAL << "Loading k-mer database " << dbFilename1 << "..." << endl;
	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile1 || !db1->deserialize(dbFile1)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename1);
	}

	LOG_NORMAL << "Loading k-mer database " << dbFilename2 << "..." << endl;
	//	auto start = std::chrono::high_resolution_clock::now();
	if (!dbFile2 || !db2->deserialize(dbFile2)) {
		throw runtime_error("Cannot open k-mer database " + dbFilename2);
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl << db1->printStats() << endl << db2->printStats() << endl;


	LOG_NORMAL << "Calculating matrix of common k-mers...";
	start = std::chrono::high_resolution_clock::now();
	SparseMatrix<uint32_t> matrix;
	calculator.db2db_sp(*db1, *db2, matrix);
	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;

	LOG_NORMAL << "Storing matrix of common k-mers in " << similarityFilename << "...";
	start = std::chrono::high_resolution_clock::now();
	ofs << "kmer-length: " << db1->getKmerLength() << " fraction: " << db1->getFraction() << " ,db-samples ,";
	std::copy(db1->getSampleNames().cbegin(), db1->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	std::copy(db2->getSampleNames().cbegin(), db2->getSampleNames().cend(), ostream_iterator<string>(ofs, ","));
	ofs << endl;

	// allocate row buffer (10000 for sample name + 100 for each row)
	char* row = new char[10000 + (db1->getSamplesCount() + db2->getSamplesCount()) * 100];
	char* ptr = row;

	ptr += sprintf(ptr, "query-samples,total-kmers,");
	ptr += num2str(db1->getSampleKmersCount().data(), db1->getSampleKmersCount().size(), ',', ptr);
	ptr += num2str(db2->getSampleKmersCount().data(), db2->getSampleKmersCount().size(), ',', ptr);
	*ptr++ = '\n';
	ofs.write(row, ptr - row);

	/*	if (sparse) {
			for(auto &row : matrix.getData())
				std::replace_if(row.begin(), row.end(),
					[below, above](pair<uint32_t, uint32_t> x) { return x.second >= below || x.second <= above; }, 0);
		}*/

	for (size_t sid = 0; sid < db1->getSamplesCount(); ++sid) {
		ptr = row;
		ptr += sprintf(ptr, "%s,%lu,", db1->getSampleNames()[sid].c_str(), db1->getSampleKmersCount()[sid]);

		/*		if (sparse) {
					ptr += matrix.saveRowSparse(sid, ptr);
				}*/
				/*		else {
							ptr += matrix.saveRow(sid, ptr);
						}*/

		*ptr++ = '\n';
		ofs.write(row, ptr - row);
	}

	size_t samples_count1 = db1->getSamplesCount();

	for (size_t sid = 0; sid < db2->getSamplesCount(); ++sid) {
		ptr = row;
		ptr += sprintf(ptr, "%s,%lu,", db2->getSampleNames()[sid].c_str(), db2->getSampleKmersCount()[sid]);

		if (params.sparseOut) {
			ptr += matrix.saveRowSparse(samples_count1 + sid, ptr);
		}
		/*		else {
					ptr += matrix.saveRow(sid, ptr);
				}*/

		*ptr++ = '\n';
		ofs.write(row, ptr - row);
	}

	delete[] row;

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;

	LOG_NORMAL << "Releasing memory...";
	start = std::chrono::high_resolution_clock::now();

	delete db1;
	delete db2;

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_NORMAL << "OK (" << dt.count() << " seconds)" << endl;
}