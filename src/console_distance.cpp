#include "console.h"

#include <chrono>
#include <cstdint>
#include <algorithm>

void DistanceConsole::run(const Params& params) {

	if (params.files.size() < 2) {
		throw usage_error(params.mode);
	}
	LOG_NORMAL("Calculating distance measures" << endl);

	const std::string& similarityFilename = params.files[0];

	std::vector<num_kmers_t> dbKmerCounts;
	uint32_t kmerLength;

	LOG_NORMAL("Loading file with common k-mer counts: " << similarityFilename << "...");
	ifstream similarityFile(similarityFilename);
	if (!similarityFile) {
		throw runtime_error("Cannot open common k-mers table: " + similarityFilename);
	}
	LOG_NORMAL("OK" << endl);

	const size_t io_buf_size = 128 << 20;
	char* io_buffer1 = new char[io_buf_size];
	char* io_buffer2 = new char[io_buf_size];
	similarityFile.rdbuf()->pubsetbuf(io_buffer1, io_buf_size);


	std::ofstream file(params.files[1]);
	std::vector<std::string> dbSampleNames;

	string tmp, in;
	double fraction;
	similarityFile >> tmp >> kmerLength >> tmp >> fraction >> tmp;

	getline(similarityFile, in); // copy sample names to output files
	if (!params.phylipOut) {
		file << "kmer-length: " << kmerLength << " fraction: " << fraction << in << endl;
	}

	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss(in);
	std::copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(), std::back_inserter(dbSampleNames));

	getline(similarityFile, in); // get number of kmers for all samples
	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss2(in);
	iss2 >> tmp >> tmp;
	std::copy(std::istream_iterator<size_t>(iss2), std::istream_iterator<size_t>(), std::back_inserter(dbKmerCounts));

	if (params.phylipOut) {
		file << dbKmerCounts.size() << endl;
	}

	std::vector<num_kmers_t> intersections_dense(dbKmerCounts.size());
	std::vector<pair<size_t, num_kmers_t>> intersections_sparse;
	std::vector<double> values_dense(dbKmerCounts.size());
	std::vector<pair<size_t, double>> values_sparse;

	const size_t bufsize = 1ULL << 30; // 1 GB  buffer
	char* outBuffer = new char[bufsize];
	char* line = new char[bufsize];
	char* begin, * end, * p;

	LOG_NORMAL("Processing rows..." << endl);
	bool triangle = false;
	bool sparseOut = params.sparseOut && !params.phylipOut; // output in Phylip format is always dense

	auto& metric = params.availableMetrics.at(params.metricName);

	// fixme: check if bufsize can be removed here - maybe some auto-adjustment of outBuffer can be applied
	for (int row_id = 0; similarityFile.getline(line, bufsize); ++row_id) {
		if ((row_id + 1) % 10 == 0) {
			LOG_NORMAL("\r" << row_id + 1 << "/" << dbKmerCounts.size() << "...                      ");
		}

		// extract name
		end = line + similarityFile.gcount() - 1; // do not count \n
		begin = line;
		p = std::find(begin, end, ',');
		string queryName(begin, p);
		begin = p + 1;

		// extract kmer count
		num_kmers_t queryKmersCount = NumericConversions::strtol(begin, &p); // assume no white characters after the number -> p points comma
		begin = p + 1;

		std::vector<num_kmers_t> queryKmerCounts(1, queryKmersCount);

		CombinedFilter<num_kmers_t> filter(
			params.metricFilters,
			params.kmerFilter,
			queryKmerCounts,
			dbKmerCounts,
			kmerLength);

		int numRead = 0;
		for (numRead = 0; end - begin > 1; ++numRead) {
			long v = NumericConversions::strtol(begin, &p);

			if (*p == ':') {
				begin = p + 1;
				num_kmers_t common = NumericConversions::strtol(begin, &p);

				if (params.phylipOut) {
					intersections_dense[v - 1] = common; // 1-based indexing in file
				}
				else {
					// sparse input always produces sparse outputs
					sparseOut = true;
					if (common > 0 && filter(common, 0, v - 1)) {
						intersections_sparse.emplace_back(v - 1, common);
					}
				}

			}
			else {
				num_kmers_t common = v;
				if (sparseOut) {
					if (common > 0 && filter(common, 0, numRead)) {
						intersections_sparse.emplace_back(numRead, common);
					}
				} else {
					// dense form
					intersections_dense[numRead] = common;
				}
			}

			begin = p + 1;
		}

		// determine if matrix is triangle
		bool emptyDiagonal = sparseOut ? (intersections_sparse.size() == 0) : (intersections_dense[0] == 0);
		if (row_id == 0 && queryName == dbSampleNames[0] && emptyDiagonal) {
			triangle = true;
		}

		// number of processed elements:
		// - dense triangle matrices - same as row id
		// - dense non-triangle - entire row
		// - others - same as input
		int numToProcess = !sparseOut
			? (triangle ? row_id : intersections_dense.size())
			: intersections_sparse.size();

		if (sparseOut) {
			values_sparse.resize(intersections_sparse.size());

			// non-empty row
			if (intersections_sparse.size() > 0) {
				
				std::transform(intersections_sparse.begin(), intersections_sparse.begin() + numToProcess, values_sparse.begin(),
					[&metric, &dbKmerCounts, queryKmersCount, kmerLength](const std::pair<size_t, num_kmers_t>& entry) {
						return std::make_pair<size_t, double>(
							entry.first + 1,
							metric(entry.second, queryKmersCount, dbKmerCounts[entry.first], kmerLength));
					});
			}
		}
		else {
			std::transform(intersections_dense.begin(), intersections_dense.begin() + numToProcess, dbKmerCounts.begin(), values_dense.begin(),
				[&metric, queryKmersCount, kmerLength](size_t intersection, num_kmers_t dbKmerCount) { 
					return metric(intersection, queryKmersCount, dbKmerCount, kmerLength); 
				});
		}
		
		char* ptr = outBuffer;
		memcpy(ptr, queryName.c_str(), queryName.size());
		ptr += queryName.size();

		if (params.phylipOut) {
			// phylip matrices are always stored in the dense form
			*ptr++ = ' ';
			ptr += num2str(values_dense.data(), numRead, ' ', ptr);
		}
		else {
			*ptr++ = ',';
			if (sparseOut) {
				for (auto& x : values_sparse) {
					ptr += num2str(x, ptr);
					*ptr++ = ',';
				}
			}
			else {
				// dense matrix - write the same number of elements as was read
				ptr += num2str(values_dense.data(), numToProcess, ',', ptr);
			}
		}

		*ptr = 0;
		size_t len = ptr - outBuffer;
		file.write(outBuffer, len);
		file << endl;
		
		if (params.phylipOut || !sparseOut) {
			intersections_dense.assign(intersections_dense.size(), 0);
		}
		else {
			intersections_sparse.clear();
		}
	}

	LOG_NORMAL("\r" << dbKmerCounts.size() << "/" << dbKmerCounts.size() << "...");
	LOG_NORMAL("OK" << endl);

	delete[] outBuffer;
	delete[] line;
	delete[] io_buffer1;
	delete[] io_buffer2;
}
