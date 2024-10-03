#include "console.h"
#include "prefix_kmer_db.h"
#include "similarity_calculator.h"

#include <chrono>
#include <cstdint>
#include <algorithm>

using sampler_t = Sampler<uint32_t, uint32_t, double>;

void All2AllPartsConsole::run(const Params& params) {
	uint32_t sampling_max_no_items = params.samplingSize;
	bool do_sampling = sampling_max_no_items != 0;
	sampler_t::strategy_t sampling_strategy = params.samplingCriterion ? sampler_t::strategy_t::best : sampler_t::strategy_t::random;

	if (params.files.size() != 2) {
		throw usage_error(params.mode);
	}

	LOG_NORMAL << "All versus all comparison (sparse computation)" << endl;
	
	const string& multipleDb = params.files[0];
	const std::string& similarityFile = params.files[1];

	SimilarityCalculator calculator(params.numThreads, params.cacheBufferMb);

	PrefixKmerDb* db_row, * db_col;

	vector<pair<string, size_t>> sample_name_count;
	vector<size_t> part_no_samples;
	uint32_t kmer_len = 0;
	double fraction = 1;

	auto t0 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t1, t2, t3, t4;

	std::chrono::duration<double> dt_load = t0 - t0;
	std::chrono::duration<double> dt_release = t0 - t0;
	std::chrono::duration<double> dt_all2all = t0 - t0;
	std::chrono::duration<double> dt_store = t0 - t0;
	std::chrono::duration<double> dt_total;

	ifstream ifs(multipleDb);

	if (!ifs) {
		throw std::runtime_error("Cannot open: " + multipleDb);
	}

	istream_iterator<string> ifs_iter(ifs);
	vector<string> input_fn(ifs_iter, istream_iterator<string>());
	uint32_t no_parts = input_fn.size();

	ifs.close();

	LOG_NORMAL << "Processing database grid of size " << no_parts << " by " << no_parts << endl;

	const size_t io_buffer_size = 32 << 20;
	char* io_buffer1 = new char[io_buffer_size];
	char* io_buffer2 = new char[io_buffer_size];
	char* io_buffer3 = new char[io_buffer_size];

	t1 = std::chrono::high_resolution_clock::now();
	for (uint32_t i = 0; i < no_parts; ++i)
	{
		db_row = new PrefixKmerDb(params.numThreads);
		string fn = input_fn[i];
		std::ifstream dbFile(fn, std::ios::binary);
		dbFile.rdbuf()->pubsetbuf(io_buffer1, io_buffer_size);

		if (!dbFile) {
			throw std::runtime_error("Cannot open k-mer database: " + fn);
		}

		db_row->deserialize(dbFile, AbstractKmerDb::DeserializationMode::SamplesOnly);
		uint32_t k = db_row->getKmerLength();
		double f = db_row->getFraction();

		if (i == 0)
		{
			kmer_len = k;
			fraction = f;
		}
		else
		{
			if (k != kmer_len) {
				throw std::runtime_error("Different k - mer lengths");
			}
			if (f != fraction) {
				throw std::runtime_error("Different fractions");
			}
		}

		auto& sample_names = db_row->getSampleNames();
		auto& kmer_counts = db_row->getSampleKmersCount();

		part_no_samples.emplace_back(sample_names.size());

		for (size_t j = 0; j < sample_names.size(); ++j)
			sample_name_count.emplace_back(sample_names[j], kmer_counts[j]);

		delete db_row;
	}
	t2 = std::chrono::high_resolution_clock::now();

	dt_load += t2 - t1;

	std::ofstream ofs(similarityFile, std::ios::binary);
	ofs.rdbuf()->pubsetbuf(io_buffer3, io_buffer_size);
	size_t line_len = 10000;

	for (const auto& x : sample_name_count)
		line_len += max<size_t>(38ull, x.first.size() + 1);			// 38 = 2*18 (max. value) + 1 (:) + 1 (,)

	char* line = new char[line_len];
	char* ptr = line;

	ofs << "kmer-length: " << kmer_len << " fraction: " << fraction << " ,db-samples ,";
	for (const auto& x : sample_name_count)
		ofs << x.first << ",";
	ofs << endl;

	ofs << "query-samples,total-kmers,";
	for (const auto& x : sample_name_count)
		ofs << x.second << ",";
	ofs << endl;

	t3 = std::chrono::high_resolution_clock::now();
	dt_store += t3 - t2;

	SparseMatrix<uint32_t>* matrix;
	vector<SparseMatrix<uint32_t>*> matrices_row;
	uint32_t k_global = 0;

	PrefixKmerDb* db_tmp = nullptr;

	sampler_t sampler(do_sampling ? sample_name_count.size() : 0, sampling_max_no_items, sampling_strategy);

	vector<uint32_t> idx_shifts = { 0 };

	for (uint32_t i_row = 0; i_row < no_parts; ++i_row)
	{
		LOG_VERBOSE << "***** Row: " << i_row + 1 << " of " << no_parts << endl;

		matrices_row.clear();
		matrices_row.resize(i_row + 1, nullptr);

		t1 = std::chrono::high_resolution_clock::now();

		db_row = new PrefixKmerDb(params.numThreads);
		string fn_row = input_fn[i_row];
		std::ifstream dbFile1(fn_row, std::ios::binary);
		dbFile1.rdbuf()->pubsetbuf(io_buffer1, io_buffer_size);
		LOG_NORMAL << "Deserializing database " << i_row + 1 << " (" << fn_row << ")" << std::flush << endl;
		db_row->deserialize(dbFile1, AbstractKmerDb::DeserializationMode::CompactedHashtables);
		uint32_t no_row_samples = db_row->getSamplesCount();

		auto cum_no_samples = idx_shifts.back();
		idx_shifts.emplace_back(cum_no_samples + no_row_samples);

		t2 = std::chrono::high_resolution_clock::now();
		dt_load += t2 - t1;

		if (i_row > 0)
		{
			uint32_t i_col = i_row - 1;
			t1 = std::chrono::high_resolution_clock::now();
			LOG_NORMAL << "Processing cell (" << i_row + 1 << "," << i_col + 1 << ")" << endl << std::flush;

			db_col = db_tmp;
			db_tmp = nullptr;

			t2 = std::chrono::high_resolution_clock::now();
			dt_load += t2 - t1;

			matrix = new SparseMatrix<uint32_t>;
			calculator.db2db_sp(*db_row, *db_col, *matrix);
			
			CombinedFilter<uint32_t> filter(
			params.metricFilters,
			params.kmerFilter,
			db_row->getSampleKmersCount(),
			db_col->getSampleKmersCount(),
			params.kmerLength);

			if (do_sampling)
			{
				matrix->add_to_sampler(filter, sampler, params.samplingCriterion, db_row->getSampleKmersCount(), db_col->getSampleKmersCount(), idx_shifts[i_row], idx_shifts[i_col], params.kmerLength, params.numThreads);
				matrix->clear();
			}
			else {
				matrix->compact(filter, params.numThreads);
			}

			t3 = std::chrono::high_resolution_clock::now();
			dt_all2all += t3 - t2;

			delete db_col;

			t4 = std::chrono::high_resolution_clock::now();
			dt_release += t4 - t3;

			matrices_row[i_col] = matrix;
		}

		for (uint32 i_col = 0; i_col + 1 < i_row; ++i_col)
		{
			t1 = std::chrono::high_resolution_clock::now();
			LOG_NORMAL << "Processing cell (" << i_row + 1 << "," << i_col + 1 << ")" << endl << std::flush;

			db_col = new PrefixKmerDb(params.numThreads);
			string fn_col = input_fn[i_col];
			std::ifstream dbFile2(fn_col, std::ios::binary);
			dbFile2.rdbuf()->pubsetbuf(io_buffer2, io_buffer_size);
			LOG_NORMAL << "Deserializing database" << i_col << " (" << fn_col << ")" << endl << std::flush;
			db_col->deserialize(dbFile2, AbstractKmerDb::DeserializationMode::CompactedHashtables);

			t2 = std::chrono::high_resolution_clock::now();
			dt_load += t2 - t1;

			matrix = new SparseMatrix<uint32_t>;
			calculator.db2db_sp(*db_row, *db_col, *matrix);

			CombinedFilter<uint32_t> filter(
				params.metricFilters,
				params.kmerFilter,
				db_row->getSampleKmersCount(),
				db_col->getSampleKmersCount(),
				params.kmerLength);

			if (do_sampling)
			{
				matrix->add_to_sampler(filter, sampler, params.samplingCriterion, db_row->getSampleKmersCount(), db_col->getSampleKmersCount(), idx_shifts[i_row], idx_shifts[i_col], params.kmerLength, params.numThreads);
				matrix->clear();
			}
			else {
				matrix->compact(filter, params.numThreads);
			}

			t3 = std::chrono::high_resolution_clock::now();
			dt_all2all += t3 - t2;

			delete db_col;

			t4 = std::chrono::high_resolution_clock::now();
			dt_release += t4 - t3;

			matrices_row[i_col] = matrix;
		}

		LOG_NORMAL << "Processing cell (" << i_row + 1 << "," << i_row + 1 << ")" << std::flush;

		t1 = std::chrono::high_resolution_clock::now();

		matrix = new SparseMatrix<uint32_t>;
		calculator.all2all_sp(*db_row, *matrix);

		CombinedFilter<uint32_t> filter(
			params.metricFilters,
			params.kmerFilter,
			db_row->getSampleKmersCount(),
			db_row->getSampleKmersCount(),
			params.kmerLength);

		if (do_sampling)
		{
			matrix->add_to_sampler(filter, sampler, params.samplingCriterion, db_row->getSampleKmersCount(), db_row->getSampleKmersCount(), idx_shifts[i_row], idx_shifts[i_row], params.kmerLength, params.numThreads);
			matrix->clear();
		}
		else {
			matrix->compact(filter, params.numThreads);
		}
		matrices_row[i_row] = matrix;

		t2 = std::chrono::high_resolution_clock::now();
		dt_all2all += t2 - t1;

		assert(db_tmp == nullptr);

		db_tmp = db_row;
		db_row = nullptr;

		t3 = std::chrono::high_resolution_clock::now();
		dt_release += t3 - t2;

		if (!do_sampling)
		{
			LOG_NORMAL << "Saving output matrix..." << endl << flush;
			for (uint32_t k = 0; k < no_row_samples; ++k)
			{
				ofs << sample_name_count[k_global + k].first << "," << sample_name_count[k_global + k].second << ",";
				ptr = line;
				size_t idx_shift = 0;

				for (uint32_t i_col = 0; i_col <= i_row; ++i_col)
				{
					matrix = matrices_row[i_col];
					ptr += matrix->saveRowSparse(k, ptr, idx_shift);
					idx_shift += part_no_samples[i_col];
				}

				ofs.write(line, ptr - line);
				ofs << endl;
			}
			LOG_NORMAL << " OK" << endl;
		}

		t4 = std::chrono::high_resolution_clock::now();
		dt_store += t4 - t3;

		k_global += no_row_samples;

		t1 = std::chrono::high_resolution_clock::now();

		for (auto& x : matrices_row)
			delete x;

		//		matrices_row.clear();

		t2 = std::chrono::high_resolution_clock::now();
		dt_release += t2 - t1;
	}

	if (do_sampling)
	{
		for (uint32_t i = 0; i < sample_name_count.size(); ++i)
		{
			ofs << sample_name_count[i].first << "," << sample_name_count[i].second << ",";

			ptr = line;
			ptr += sampler.saveRowSparse(i, ptr, 0);

			ofs.write(line, ptr - line);
			ofs << endl;
		}
	}

	t1 = std::chrono::high_resolution_clock::now();
	ofs.close();
	t2 = std::chrono::high_resolution_clock::now();
	dt_store += t2 - t1;

	LOG_NORMAL << "Database grid procesed successfully" << endl;
	LOG_NORMAL << "  Load time   : " << dt_load.count() << " seconds" << endl;
	LOG_NORMAL << "  All2All time: " << dt_all2all.count() << " seconds" << endl;
	LOG_NORMAL << "  Store time  : " << dt_store.count() << " seconds" << endl;

	LOG_NORMAL << "Releasing memory...";
	t2 = std::chrono::high_resolution_clock::now();
	delete[] line;
	delete[] io_buffer1;
	delete[] io_buffer2;
	delete[] io_buffer3;
	delete db_tmp;
	t3 = std::chrono::high_resolution_clock::now();
	dt_release = t3 - t2;
	LOG_NORMAL << "OK (" << dt_release.count() << " seconds)" << endl;

	dt_total = t3 - t0;
	LOG_NORMAL << "Total time  : " << dt_total.count() << " seconds" << endl;
}