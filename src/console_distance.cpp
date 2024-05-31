#include "console.h"

#include <chrono>
#include <cstdint>
#include <algorithm>

void DistanceConsole::run(const Params& params) {

	if (params.files.size() < 1) {
		throw usage_error(params.mode);
	}
	LOG_NORMAL << "Calculating distance measures" << endl;

	const std::string& similarityFilename = params.files[0];

	std::vector<size_t> kmersCount;
	uint32_t kmerLength;

	LOG_NORMAL << "Loading file with common k-mer counts: " << similarityFilename << "...";
	ifstream similarityFile(similarityFilename);
	if (!similarityFile) {
		throw runtime_error("Cannot open common k-mers table: " + similarityFilename);
	}
	LOG_NORMAL << "OK" << endl;

	const size_t io_buf_size = 128 << 20;
	char* io_buffer1 = new char[io_buf_size];
	char* io_buffer2 = new char[io_buf_size];
	similarityFile.rdbuf()->pubsetbuf(io_buffer1, io_buf_size);


	std::vector<std::ofstream> files(params.metrics.size());
	std::vector<std::string> dbSampleNames;

	for (size_t i = 0; i < files.size(); ++i) {
		files[i].open(similarityFilename + "." + params.metrics[i].first);
	}

	string tmp, in;
	double fraction;
	similarityFile >> tmp >> kmerLength >> tmp >> fraction >> tmp;

	getline(similarityFile, in); // copy sample names to output files
	if (!params.phylipOut) {
		for (auto& f : files) {
			f << "kmer-length: " << kmerLength << " fraction: " << fraction << in << endl;
		}
	}

	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss(in);
	std::copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(), std::back_inserter(dbSampleNames));

	getline(similarityFile, in); // get number of kmers for all samples
	std::replace(in.begin(), in.end(), ',', ' ');
	istringstream iss2(in);
	iss2 >> tmp >> tmp;
	std::copy(std::istream_iterator<size_t>(iss2), std::istream_iterator<size_t>(), std::back_inserter(kmersCount));

	if (params.phylipOut) {
		for (auto& f : files) {
			f << kmersCount.size() << endl;
		}
	}

	std::vector<size_t> intersections_dense(kmersCount.size());
	std::vector<pair<size_t, size_t>> intersections_sparse;
	std::vector<double> values_dense(kmersCount.size());
	std::vector<pair<size_t, double>> values_sparse;

	const size_t bufsize = 1ULL << 30; // 1 GB  buffer
	char* outBuffer = new char[bufsize];
	char* line = new char[bufsize];
	char* begin, * end, * p;

	LOG_NORMAL << "Processing rows..." << endl;
	bool triangle = false;
	bool sparseIn = false, sparseOut = params.sparseOut;

	// fixme: check if bufsize can be removed here - maybe some auto-adjustment of outBuffer can be applied
	for (int row_id = 0; similarityFile.getline(line, bufsize); ++row_id) {
		if ((row_id + 1) % 10 == 0) {
			LOG_NORMAL << "\r" << row_id + 1 << "/" << kmersCount.size() << "...                      " << std::flush;
		}

		// extract name
		end = line + similarityFile.gcount() - 1; // do not count \n
		begin = line;
		p = std::find(begin, end, ',');
		string queryName(begin, p);
		begin = p + 1;

		// extract kmer count
		uint64_t queryKmersCount = NumericConversions::strtol(begin, &p); // assume no white characters after the number -> p points comma
		begin = p + 1;

		int numRead = 0;
		for (numRead = 0; end - begin > 1; ++numRead) {
			long v = NumericConversions::strtol(begin, &p);

			if (*p == ':') {
				// sparse input always produces sparse outputs
				sparseIn = sparseOut = true;

				begin = p + 1;
				long common = NumericConversions::strtol(begin, &p);

				if (params.phylipOut || !sparseOut)
					intersections_dense[v - 1] = common; // 1-based indexing in file
				else
					intersections_sparse.emplace_back(v - 1, common);

			}
			else {
				// dense form
				intersections_dense[numRead] = v;
			}

			begin = p + 1;
		}

		// determine if matrix is triangle
		if (row_id == 0 && queryName == dbSampleNames[0] && intersections_dense[0] == 0) {
			triangle = true;
		}

		for (size_t m = 0; m < params.metrics.size(); ++m) {
			auto& metric = params.metrics[m].second;

			// number of processed elements:
			// - triangle matrices - same as row id
			// - non-triangle sparse matrices - entire row 
			// - others - same as input
			int numToProcess = triangle ? row_id : (sparseIn ? intersections_dense.size() : numRead);

			if (sparseIn && !params.phylipOut && sparseOut)
			{
				values_sparse.clear();
				for (size_t i = 0; i < intersections_sparse.size(); ++i)
				{
					size_t id = intersections_sparse[i].first;
					values_sparse.emplace_back(id + 1, metric(intersections_sparse[i].second, queryKmersCount, kmersCount[id], kmerLength));
				}
			}
			else
				std::transform(intersections_dense.begin(), intersections_dense.begin() + numToProcess, kmersCount.begin(), values_dense.begin(),
					[&metric, queryKmersCount, kmerLength](size_t intersection, size_t dbKmerCount)->double { return  metric(intersection, queryKmersCount, dbKmerCount, kmerLength); });

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
					if (sparseIn)
					{
						for (auto& x : values_sparse)
							if (x.second < params.below && x.second > params.above)
							{
								ptr += num2str(x, ptr);
								*ptr++ = ',';
							}
					}
					else
					{
						// sparse value
						double sparse_val = 0.0;

						// if boundaries are specified - remove elements outside boundaries
						if (params.below != params.MAX_BELOW || params.above != params.MIN_ABOVE) {
							sparse_val = (params.below != params.MAX_BELOW) ? params.MAX_BELOW : params.MIN_ABOVE;

							std::replace_if(values_dense.begin(), values_dense.end(),
								[&params](double x) { return x >= params.below || x <= params.above; }, sparse_val);
						}

						ptr += num2str_sparse(values_dense.data(), numToProcess, ',', ptr, sparse_val);
					}
				}
				else {
					// dense matrix - write the same number of elements as was read
					ptr += num2str(values_dense.data(), numToProcess, ',', ptr);
				}
			}

			*ptr = 0;
			size_t len = ptr - outBuffer;
			files[m].write(outBuffer, len);
			files[m] << endl;
		}

		if (sparseIn) {
			if (params.phylipOut || !sparseOut)
				intersections_dense.assign(intersections_dense.size(), 0);
			else
				intersections_sparse.clear();
		}
	}

	LOG_NORMAL << "\r" << kmersCount.size() << "/" << kmersCount.size() << "...";
	LOG_NORMAL << "OK" << endl;

	delete[] outBuffer;
	delete[] line;
	delete[] io_buffer1;
	delete[] io_buffer2;
}
