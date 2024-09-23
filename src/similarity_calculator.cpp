// build -k 25 -t 32 -multisample-fasta fl.txt part_00000.kdb
// all2all-sp -sparse -t 32 part_00000.kdb part_00000.a2as

#include "similarity_calculator.h"
#include "parallel_sorter.h"
#include "../libs/refresh/sort/lib/pdqsort_par.h"

#if defined(ARCH_X64)
#include <emmintrin.h>
#include <immintrin.h>
#elif defined(ARCH_ARM)
#include <arm_neon.h>
#endif

#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "instr_set_detect.h"
#include <barrier>
#include <future>

// *****************************************************************************************
//
SimilarityCalculator::SimilarityCalculator(int _num_threads, size_t cacheBufferMb) : 
	// hardware_concurrency may return 0
	num_threads(_num_threads > 0 ? _num_threads : std::max((int)std::thread::hardware_concurrency(), 1)),
	cacheBufferMb(cacheBufferMb),
	atp(4, 1024)
{

#if defined(ARCH_X64)
	// !!! TODO: do zmiany ten instrset_detect()
	auto instr = InstrSetDetect::GetInstr();
	avx2_present = (instr == InstrSetDetect::Instr::AVX2);
#elif defined(ARCH_ARM)
	// 
#endif
}

// *****************************************************************************************
//
void SimilarityCalculator::all2all(PrefixKmerDb& db, LowerTriangularMatrix<uint32_t>& matrix) const 
{
	// get stuff from database
	auto& patterns = db.getPatterns();
	int samples_count = db.getSamplesCount();

	matrix.resize(samples_count);
	//matrix.clear();

	size_t bufsize = cacheBufferMb * 1000000 / sizeof(uint32_t);
	std::vector<uint32_t> patternsBuffer(bufsize);

	std::vector<uint32_t*> rawPatterns(bufsize);
	std::vector<std::tuple<sample_id_t, uint32_t, uint32_t>> sample2pattern(bufsize);

	std::vector<std::thread> workers_matrix(num_threads);
	std::vector<std::thread> workers_decomp(num_threads);
	std::vector<std::thread> workers_sample2patterns(num_threads);
	std::vector<std::thread> workers_histogram(num_threads);

	int first_pid;

	// Set number of k-mers to internal nodes
	int num_patterns = patterns.size();
	for (int i = num_patterns - 1; i > 0; --i)
	{
		int parent_id = patterns[i].get_parent_id();

		if (parent_id >= 0)
			patterns[parent_id].add_num_kmers(patterns[i].get_num_kmers());
	}

	// determine ranges of blocks processed by threads 
	int no_ranges = num_threads * 16;
	std::vector<size_t> workerRanges(no_ranges + 1, sample2pattern.size());
	std::vector<int> v_range_ids;

	for (int i = no_ranges - 1; i >= 0; --i)
		v_range_ids.push_back(i);

	RegisteringQueue<int> tasks_matrix_queue(1);
	Semaphore semaphore_matrix;
	std::atomic<uint64_t> numAdditions(0);

	// ranges of patterns to decompress
	std::vector<int> v_range_patterns(no_ranges + 1);
	RegisteringQueue<pair<int, uint32_t>> tasks_decomp_queue(1);
	Semaphore semaphore_decomp;

	RegisteringQueue<int> tasks_sample2patterns_queue(1);
	Semaphore semaphore_sample2patterns;

	Semaphore semaphore_hist_first;		// semaphore for first stage in histogram calculation
	Semaphore semaphore_hist_second;		// semaphore for second stage in histogram calculation

	int no_hist_parts = num_threads * 1;
	std::vector<std::vector<int>> hist_sample_ids(no_hist_parts, std::vector<int>(samples_count, 0));
	RegisteringQueue<int> tasks_histogram_queue(1);
	std::vector<int> hist_boundary_values(num_threads, 0);
	std::vector<pair<int, uint32_t>> v_tmp;
	std::vector<int> v_tmp_int;
	std::vector<int> v_hist_ids(num_threads);
	for (int i = 0; i < num_threads; ++i)
		v_hist_ids[i] = i;

	bool use_busy_wait = num_threads >= 8;

	// Decompress patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_decomp[tid] = std::thread([this, samples_count, &patterns, &sample2pattern, &v_range_patterns, &tasks_decomp_queue, &semaphore_decomp, &rawPatterns, &first_pid, &patternsBuffer, &hist_sample_ids]
		{
			pair<int, uint32_t> decomp_task;

			while (tasks_decomp_queue.Pop(decomp_task))
			{
				int part_id = decomp_task.first;
				auto &my_hist_sample_ids = hist_sample_ids[part_id];
				my_hist_sample_ids.clear();
				my_hist_sample_ids.resize(samples_count, 0);

				int f_pid = v_range_patterns[part_id];
				int l_pid = v_range_patterns[part_id + 1];
				auto currentPtr = patternsBuffer.data() + decomp_task.second;

				for (int pid = f_pid; pid < l_pid; ++pid)
				{
					const auto& pattern = patterns[pid];

					currentPtr += pattern.get_num_samples();
					uint32_t* out = currentPtr;		// start from the end

													// decode all samples from pattern and its parents
					int64_t current_id = pid;
					while (current_id >= 0) {
						const auto& cur = patterns[current_id];
						auto parent_id = cur.get_parent_id();

						if (parent_id >= 0)
//							_mm_prefetch((const char*)(patterns.data() + parent_id), _MM_HINT_T0);
#ifdef WIN32
							_mm_prefetch((const char*)(patterns.data() + parent_id), _MM_HINT_T0);
#else
							__builtin_prefetch(patterns.data() + parent_id);
#endif

						out -= cur.get_num_local_samples();
						cur.decodeSamples(out);

						current_id = parent_id;
					}
					rawPatterns[pid - first_pid] = out; // begin of unpacked pattern

					int num_samples = pattern.get_num_samples();
					int num_local_samples = pattern.get_num_local_samples();
					for (int i = num_samples - num_local_samples; i < num_samples; ++i)
						++my_hist_sample_ids[out[i]];
				}

				semaphore_decomp.dec();
			}
		});
	}

	// Sample2patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_sample2patterns[tid] = std::thread([
			this, &patterns, &sample2pattern, &v_range_patterns, &tasks_sample2patterns_queue, 
				&semaphore_sample2patterns, &rawPatterns, &first_pid, &patternsBuffer, &hist_sample_ids]
		{
//			while (!tasks_sample2patterns_queue.IsCompleted())
//			{
				int part_id;

//				if (tasks_sample2patterns_queue.Pop(part_id))
				while (tasks_sample2patterns_queue.Pop(part_id))
				{
					auto &my_hist_sample_ids = hist_sample_ids[part_id];

					int f_pid = v_range_patterns[part_id];
					int l_pid = v_range_patterns[part_id + 1];

					for (int pid = f_pid; pid < l_pid; ++pid)
					{
						const auto& pattern = patterns[pid];

						auto out = rawPatterns[pid - first_pid];
						int num_samples = pattern.get_num_samples();
						int num_local_samples = pattern.get_num_local_samples();

						for (int i = num_samples - num_local_samples; i < num_samples; ++i)
						{
							int pos = my_hist_sample_ids[out[i]]++;

							sample2pattern[pos] = make_tuple(out[i], pid - first_pid, i);
						}
					}

					semaphore_sample2patterns.dec();
				}
//			}
		});
	}

	// increment array elements in threads
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_matrix[tid] = std::thread([this, &patterns, &sample2pattern, &rawPatterns, &workerRanges, &matrix, &tasks_matrix_queue, &first_pid, &semaphore_matrix, &numAdditions] {
			
			int range_id;

			while (tasks_matrix_queue.Pop(range_id))
			{
				for (size_t id = workerRanges[range_id]; id < workerRanges[range_id + 1]; ) {
					sample_id_t Si = get<0>(sample2pattern[id]);
					uint32_t *row = matrix[Si];
					while (id < workerRanges[range_id + 1] && get<0>(sample2pattern[id]) == Si) {
						int local_pid = get<1>(sample2pattern[id]);
						const auto& pattern = patterns[local_pid + first_pid];

						uint32_t* rawData = rawPatterns[local_pid];
						int num_samples = get<2>(sample2pattern[id]);
						uint32_t to_add = pattern.get_num_kmers();

#ifdef ALL_STATS
						localAdditions += num_samples;
#endif

						auto *p = rawData;

						row_add(row, p, num_samples, to_add, avx2_present);
						++id;
					}
				}
				semaphore_matrix.dec();
			}

#ifdef ALL_STATS
			numAdditions.fetch_add(localAdditions);
#endif
		});
	}

	// calculate histogram
	for (int tid = 0; tid < num_threads; ++tid)
		workers_histogram[tid] = std::thread([&hist_boundary_values, &tasks_histogram_queue, samples_count, this, &v_tmp,
			&hist_sample_ids, &semaphore_hist_first, &semaphore_hist_second, &use_busy_wait] {
		int task_id;

		while (tasks_histogram_queue.Pop(task_id))
		{
			int first_sample = samples_count * task_id / num_threads;
			int last_sample = samples_count * (task_id + 1) / num_threads;

			int sum = 0;
			int max_j = v_tmp.size();
			for (int i = first_sample; i < last_sample; ++i)
				for (int j = 0; j < max_j; ++j)
				{
					int x = hist_sample_ids[j][i];
					hist_sample_ids[j][i] = sum;
					sum += x;
				}

			hist_boundary_values[task_id] = sum;
			semaphore_hist_first.dec_notify_all();
			if(use_busy_wait)
				semaphore_hist_first.waitForZeroBusy();
			else
				semaphore_hist_first.waitForZero();

			int to_add = 0;
			for (int i = 0; i < task_id; ++i)
				to_add += hist_boundary_values[i];

			for (int j = 0; j < max_j; ++j)
				for (int i = first_sample; i < last_sample; ++i)
					hist_sample_ids[j][i] += to_add;
			semaphore_hist_second.dec();
		}
	});

	// process all patterns in blocks determined by buffer size
	//LOG_NORMAL << std::endl;

	double decomp_time = 0;
	double hist_time = 0;
	double patterns_time = 0;

	size_t percent = 0;
	for (size_t pid = 0; pid < patterns.size(); ) {
		
		if (pid * 100 / patterns.size() >= percent) {
			LOG_NORMAL << "\r" << percent << "%" << std::flush;
			++percent;
		}

		first_pid = pid;
		size_t samplesCount = 0;

		auto t1 = std::chrono::high_resolution_clock::now();

		size_t part_size = bufsize / no_hist_parts + 1;
		size_t next_boundary = 0;
		int part_id = 0;

		v_tmp.clear();
		v_tmp_int.clear();

		// !!! Fix me: what if no_samples > bufsize?
		while (pid < patterns.size() && samplesCount + patterns[pid].get_num_samples() < bufsize && part_id < no_hist_parts) {
			const auto& pattern = patterns[pid];

			if (samplesCount >= next_boundary)
			{
				next_boundary += part_size;
				v_range_patterns[part_id] = pid;
				v_tmp.push_back(make_pair(part_id, (uint32_t)samplesCount));
				v_tmp_int.push_back(part_id);

				++part_id;
			}

			samplesCount += pattern.get_num_samples();
			++pid;
		}

		v_range_patterns[part_id] = pid;
		semaphore_decomp.inc(v_tmp.size());
		tasks_decomp_queue.PushRange(v_tmp);

		if(use_busy_wait)
			semaphore_decomp.waitForZeroBusy();
		else
			semaphore_decomp.waitForZero();

		auto t2 = std::chrono::high_resolution_clock::now();

		semaphore_hist_first.inc(num_threads);
		semaphore_hist_second.inc(num_threads);
		tasks_histogram_queue.PushRange(v_hist_ids);

		if (use_busy_wait)
			semaphore_hist_second.waitForZeroBusy();
		else
			semaphore_hist_second.waitForZero();

		auto t3 = std::chrono::high_resolution_clock::now();

		// generate sample to pattern mapping
		int num_items = std::accumulate(hist_boundary_values.begin(), hist_boundary_values.end(), 0);
		sample2pattern.resize(num_items);		// resize to number of items in the current fragment

		semaphore_sample2patterns.inc(v_tmp_int.size());
		tasks_sample2patterns_queue.PushRange(v_tmp_int);

		if (use_busy_wait)
			semaphore_sample2patterns.waitForZeroBusy();
		else
			semaphore_sample2patterns.waitForZero();

		auto t4 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> dt1 = t2 - t1;
		std::chrono::duration<double> dt2 = t3 - t2;
		std::chrono::duration<double> dt3 = t4 - t3;

		decomp_time += dt1.count();
		hist_time += dt2.count();
		patterns_time += dt3.count();
		
		// determine ranges of blocks processed by threads 
		workerRanges.assign(no_ranges + 1, sample2pattern.size());
		workerRanges[0] = 0;
		size_t workerBlock = sample2pattern.size() / no_ranges;

		auto currentIndex = workerBlock;

		if (currentIndex > 0) {
			for (int rid = 0; rid < no_ranges - 1; ++rid) {
				// make sure no single row of matrix is updated by multiple workers
				auto it = std::upper_bound(
					sample2pattern.begin() + currentIndex,
					sample2pattern.end(),
					*(sample2pattern.begin() + currentIndex - 1),
					[](const std::tuple<sample_id_t, uint32_t, uint32_t>& x, const std::tuple<sample_id_t, uint32_t, uint32_t>& y) -> bool {
					return get<0>(x) < get<0>(y); });	// Necessary, because we are looking for end of sample data

				size_t range = it - sample2pattern.begin();
				workerRanges[rid + 1] = range;
				currentIndex = range + workerBlock;

				if (currentIndex >= sample2pattern.size())
					break;
			}
		}

		// this should never happen
		if (workerRanges[no_ranges] != sample2pattern.size()) {
			throw std::runtime_error("ERROR in FastKmerDb::calculateSimilarity(): Invalid ranges");
		}
		

		semaphore_matrix.inc(no_ranges);
		tasks_matrix_queue.PushRange(v_range_ids);

		if (use_busy_wait)
			semaphore_matrix.waitForZeroBusy();
		else
			semaphore_matrix.waitForZero();
	}
	
	tasks_matrix_queue.MarkCompleted();
	tasks_decomp_queue.MarkCompleted();
	tasks_sample2patterns_queue.MarkCompleted();
	tasks_histogram_queue.MarkCompleted();

	for (auto & w : workers_matrix)
		w.join();

	for (auto & w : workers_decomp)
		w.join();

	for (auto & w : workers_sample2patterns)
		w.join();

	for (auto &w : workers_histogram)
		w.join();

	LOG_NORMAL << "\r100%" << std::flush << "\r" << std::flush;

	LOG_VERBOSE << "  Decomp   time: " << decomp_time << endl;
	LOG_VERBOSE << "  Hist     time: " << hist_time << endl;
	LOG_VERBOSE << "  Sample2p time: " << patterns_time << endl;

#ifdef ALL_STATS
	LOG_NORMAL << "Number of additions:" << numAdditions << endl;
#endif
}

// *****************************************************************************************
//
void SimilarityCalculator::all2all_sp(PrefixKmerDb& db, SparseMatrix<uint32_t>& matrix) const
{
	auto& patterns = db.getPatterns();
	
	matrix.resize(db.getSamplesCount());

	size_t buf_size = cacheBufferMb * 1000000 / sizeof(uint32_t);

	if (buf_size < db.getSamplesCount())
		buf_size = 2 * db.getSamplesCount();

	std::vector<uint32_t> patternsBuffer_cur(buf_size);
	std::vector<uint32_t> patternsBuffer_next(buf_size);

	// Set number of k-mers to internal nodes
//	int num_patterns = patterns.size();
/*	for (int i = num_patterns - 1; i > 0; --i)
	{
		int parent_id = patterns[i].get_parent_id();

		if (parent_id >= 0)
			patterns[parent_id].add_num_kmers(patterns[i].get_num_kmers());
	}*/

	vector<pair<uint32_t, uint32_t*>> pattern_ptrs_cur;
	vector<pair<uint32_t, uint32_t*>> pattern_ptrs_next;

	auto t1 = std::chrono::high_resolution_clock::now();

	barrier bar_patterns(num_threads + 1);

	size_t percent = 0;

	thread thr_pattern_decompress([&] {
		auto& patterns = db.getPatterns();
		size_t buf_pos = 0;

	
		for (size_t i = 0; i < patterns.size(); ++i)
		{
			if (i * 100 / patterns.size() >= percent) {
				LOG_NORMAL << "\r" << percent << "%" << std::flush;
				++percent;
			}
			
			const auto& ps = patterns[i];

			if (buf_pos + ps.get_num_samples() >= buf_size)
			{
				bar_patterns.arrive_and_wait();
				swap(pattern_ptrs_cur, pattern_ptrs_next);
				swap(patternsBuffer_cur, patternsBuffer_next);
				bar_patterns.arrive_and_wait();

				pattern_ptrs_next.clear();
				buf_pos = 0;
			}

			uint32_t* buf_ptr = patternsBuffer_next.data() + buf_pos;
			decode_pattern_samples(patterns, i, buf_ptr);

			pattern_ptrs_next.emplace_back(ps.get_num_samples(), buf_ptr);
			buf_pos += ps.get_num_samples();
		}

		bar_patterns.arrive_and_wait();
		swap(pattern_ptrs_cur, pattern_ptrs_next);
		swap(patternsBuffer_cur, patternsBuffer_next);
		bar_patterns.arrive_and_wait();

		if (!patternsBuffer_cur.empty())
		{
			bar_patterns.arrive_and_wait();
			patternsBuffer_cur.clear();
			bar_patterns.arrive_and_wait();
		}

		});

	vector<thread> thr_workers;
	thr_workers.reserve(num_threads);

	for (int i = 0; i < num_threads; ++i)
		thr_workers.emplace_back([&, i] {
		uint32_t thread_id = i;
		uint32_t my_mod1 = thread_id % num_threads;
		uint32_t my_mod2 = 2 * num_threads - 1 - my_mod1;
		uint32_t pid = 0;

		bar_patterns.arrive_and_wait();
		bar_patterns.arrive_and_wait();

		vector<uint32_t> my_js;
		vector<pair<uint32_t, uint32_t>> hash_j;

		const uint32_t pf_skip = 32;

		while (true)
		{
			if (patternsBuffer_cur.empty())
				break;

//			LOG_NORMAL << to_string(thread_id) + " ";

			for (const auto& pat : pattern_ptrs_cur)
			{
				uint32_t to_add = patterns[pid++].get_num_kmers();

//				if(thread_id == 0)
//					cout << pat.first << " ";

#if 0	// OLD
				for (uint32_t j = 0; j < pat.first; ++j)
				{
					uint32_t* pat_data = pat.second;
					uint32_t row_id = pat_data[j];

					uint32_t row_id_nt = row_id % (2 * num_threads);
					if (row_id_nt != my_mod1 && row_id_nt != my_mod2)
						continue;

					auto& row_data = matrix[row_id];

					for (uint32_t k = 0; k < j; ++k)
						row_data[pat_data[k]] += to_add;
				}

#else
				my_js.clear();

				uint32_t* pat_data = pat.second;

				for (uint32_t j = 0; j < pat.first; ++j)
				{
					uint32_t row_id_nt = pat_data[j] % (2 * num_threads);
					if (row_id_nt == my_mod1 || row_id_nt == my_mod2)
						my_js.emplace_back(j);
				}

				for (size_t jj = 0; jj < my_js.size(); ++jj)
				{
					if (my_js[jj] < pf_skip)
					{
						if (jj + 1 < my_js.size())
						{
							uint32_t j_next = my_js[jj + 1];
							uint32_t row_id_next = pat_data[j_next];
							auto& row_data_next = matrix[row_id_next];

							for (uint32_t k = 0; k < j_next; ++k)
								row_data_next.prefetch(pat_data[k]);
						}

						uint32_t j = my_js[jj];
						uint32_t row_id = pat_data[j];

						auto& row_data = matrix[row_id];

						for (uint32_t k = 0; k < j; ++k)
							row_data[pat_data[k]] += to_add;
					}
					else
					{
						uint32_t j = my_js[jj];
						uint32_t row_id = pat_data[j];
						uint32_t row_id_next = (jj + 1 < my_js.size()) ? pat_data[my_js[jj + 1]] : 0;

						auto& row_data = matrix[row_id];

						for (uint32_t k = 0; k < j - pf_skip; ++k)
						{
							row_data.prefetch(pat_data[k+pf_skip]);
							row_data[pat_data[k]] += to_add;
						}

						for (uint32_t k = j - pf_skip; k < j; ++k)
						{
							matrix[row_id_next].prefetch(pat_data[k - (j - pf_skip)]);			// for simplicity (to avoid if) can prefetch row_id_next=0
							row_data[pat_data[k]] += to_add;
						}
					}
				}
#endif
			}

			bar_patterns.arrive_and_wait();
			bar_patterns.arrive_and_wait();
		}
			});

	thr_pattern_decompress.join();
	for (auto& t : thr_workers)
		t.join();

	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dt1 = t2 - t1;

	LOG_NORMAL << "\r100%" << std::flush << "\r" << std::flush;

	//LOG_VERBOSE << "All2All Sparse computation time: " << dt1.count() << endl;
}

// *****************************************************************************************
//
template <>
void  SimilarityCalculator::one2all<true>(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<uint32_t>& similarities) const {
	
	// get stuff from database
	const auto& patterns = db.getPatterns();
	int samples_count = db.getSamplesCount();
	const auto& hashtables = db.getHashtables();
	
	similarities.resize(samples_count, 0);

	std::vector<std::vector<uint32_t>> localSimilarities(num_threads, std::vector<uint32_t>(samples_count));

	std::unordered_map<pattern_id_t, int32_t> patterns2count;

	std::chrono::duration<double> dt;
	auto start = std::chrono::high_resolution_clock::now();

	// iterate over kmers in analyzed sample
	for (size_t i = 0; i < kmersCount; ++i) {

		if (i + PREFETCH_DIST < kmersCount) {
			kmer_t prefetch_kmer = kmers[i + PREFETCH_DIST];

			kmer_t prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
			suffix_t suffix = GET_SUFFIX(prefetch_kmer);

			hashtables[prefix].prefetch(suffix);
		}

		// check if kmer exists in a database
		kmer_t kmer = kmers[i];
		kmer_t prefix = GET_PREFIX_SHIFTED(kmer);
		suffix_t suffix = GET_SUFFIX(kmer);

		auto entry = hashtables[prefix].find(suffix);

		if (entry != nullptr) {
			auto pid = *entry;
			const auto& pattern = patterns[pid];

			if (pattern.get_num_kmers() == 0)
				continue;

			++patterns2count[pid];
		}
	}

	std::vector<std::pair<pattern_id_t, int32_t>> patterns2countVector(patterns2count.size());
	int i = 0;
	uint64_t tot_pattern_lengths = 0;
	for (const auto& entry : patterns2count) {
		tot_pattern_lengths += patterns[entry.first].get_num_samples();
	}

	vector<uint32_t> range_boundaries(num_threads + 1, patterns2count.size());
	uint64_t sum_pattern_lengths = 0;
	uint64_t part_length = (tot_pattern_lengths + num_threads - 1) / num_threads;
	uint32_t range_id = 1;
	for (auto &entry : patterns2count)
	{
		sum_pattern_lengths += patterns[entry.first].get_num_samples();
		patterns2countVector[i++] = entry;

		if (sum_pattern_lengths > range_id * part_length)
			range_boundaries[range_id++] = i;
	}
	range_boundaries[0] = 0;

	patterns2count.clear();


	dt = std::chrono::high_resolution_clock::now() - start;
	//LOG_VERBOSE << "Pattern listing time: " << dt.count() << endl ;

	start = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> workers(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, &patterns, samples_count, tid, &patterns2countVector, &localSimilarities, &range_boundaries]() {
			std::vector<uint32_t> samples(samples_count);

			size_t lo = range_boundaries[tid];
			size_t hi = range_boundaries[tid + 1];
			auto &my_localSimilarities = localSimilarities[tid];

			for (size_t id = lo; id < hi; ++id) {
				if (id + 1 < hi)
//					_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);
#ifdef WIN32
					_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);
#else
					__builtin_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first));
#endif

				int to_add = patterns2countVector[id].second;
				auto pid = patterns2countVector[id].first;

				int num_samples = decode_pattern_samples(patterns, pid, samples.data());
/*				const auto& pattern = patterns[pid];
				int num_samples = pattern.get_num_samples();

				uint32_t* out = samples.data() + pattern.get_num_samples(); // start from the end

				int64_t current_id = pid;
				while (current_id >= 0) {
					const auto& cur = patterns[current_id];

					out -= cur.get_num_local_samples();
					cur.decodeSamples(out);

					current_id = cur.get_parent_id();
				}*/

				auto *p = samples.data();

				int i;
				for (i = 0; i + 4 <= num_samples; i += 4)
				{
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
					my_localSimilarities[*p++] += to_add;
				}
				num_samples -= i;

				switch (num_samples)
				{
				case 3: my_localSimilarities[*p++] += to_add;
				case 2: my_localSimilarities[*p++] += to_add;
				case 1: my_localSimilarities[*p++] += to_add;
				}

			}
		});
	}

	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid].join();
		std::transform(localSimilarities[tid].begin(), localSimilarities[tid].end(), similarities.begin(), similarities.begin(), [](uint32_t a, uint32_t b)->uint32_t {
			return a + b;
		});
	}

	dt = std::chrono::high_resolution_clock::now() - start;
	//LOG_VERBOSE << "Pattern unpacking time: " << dt.count() << endl ;
}

// *****************************************************************************************
//
template <>
void SimilarityCalculator::one2all<false>(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<uint32_t>& similarities) const
{
	// get stuff from database
	const auto& patterns = db.getPatterns();
	int samples_count = db.getSamplesCount();
	const auto& hashtables = db.getHashtables();

	similarities.resize(samples_count, 0);

	std::unordered_map<pattern_id_t, int32_t> patterns2count;

	std::chrono::duration<double> dt;
	auto start = std::chrono::high_resolution_clock::now();

	// iterate over kmers in analyzed sample
	for (size_t i = 0; i < kmersCount; ++i) {

		if (i + PREFETCH_DIST < kmersCount) {
			kmer_t prefetch_kmer = kmers[i + PREFETCH_DIST];

			kmer_t prefix = GET_PREFIX_SHIFTED(prefetch_kmer);
			suffix_t suffix = GET_SUFFIX(prefetch_kmer);

			hashtables[prefix].prefetch(suffix);
		}

		// check if kmer exists in a database
		kmer_t kmer = kmers[i];
		kmer_t prefix = GET_PREFIX_SHIFTED(kmer);
		suffix_t suffix = GET_SUFFIX(kmer);

		auto entry = hashtables[prefix].find(suffix);

		if (entry != nullptr) {
			auto pid = *entry;
			const auto& pattern = patterns[pid];

			if (pattern.get_num_kmers() == 0)
				continue;

			++patterns2count[pid];
		}
	}

	std::vector<std::pair<pattern_id_t, int32_t>> patterns2countVector(patterns2count.size());
	int i = 0;
	
	uint64_t sum_pattern_lengths = 0;
	for (auto &entry : patterns2count)
	{
		sum_pattern_lengths += patterns[entry.first].get_num_samples();
		patterns2countVector[i++] = entry;	
	}
	
	patterns2count.clear();

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_VERBOSE << "Pattern listing time: " << dt.count() << endl ;

	start = std::chrono::high_resolution_clock::now();
	
	std::vector<uint32_t> samples(samples_count);

	size_t lo = 0;
	size_t hi = patterns2countVector.size();
		
	for (size_t id = lo; id < hi; ++id) {
		if (id + 1 < hi)
//			_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);
#ifdef WIN32
			_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);
#else
			__builtin_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first));
#endif

		auto pid = patterns2countVector[id].first;
		const auto& pattern = patterns[pid];
		int num_samples = pattern.get_num_samples();
		int to_add = patterns2countVector[id].second;

		uint32_t* out = samples.data() + pattern.get_num_samples(); // start from the end

		int64_t current_id = pid;
		while (current_id >= 0) {
			const auto& cur = patterns[current_id];

			out -= cur.get_num_local_samples();
			cur.decodeSamples(out);

			current_id = cur.get_parent_id();
		}

		auto *p = samples.data();

		int i;
		for (i = 0; i + 4 <= num_samples; i += 4)
		{
			similarities[*p++] += to_add;
			similarities[*p++] += to_add;
			similarities[*p++] += to_add;
			similarities[*p++] += to_add;
		}
		num_samples -= i;

		switch (num_samples)
		{
		case 3: similarities[*p++] += to_add;
		case 2: similarities[*p++] += to_add;
		case 1: similarities[*p++] += to_add;
		}

	}

	dt = std::chrono::high_resolution_clock::now() - start;
	LOG_VERBOSE << "Pattern unpacking time: " << dt.count() << endl ;
}

// *****************************************************************************************
//
void SimilarityCalculator::db2db(const PrefixKmerDb& db1, const PrefixKmerDb& db2, LowerTriangularMatrix<uint32_t>& matrix) const
{
	// Just part of the code

/*	auto& ht1 = db1.getHashtables();
	auto& ht2 = db2.getHashtables();

	vector<thread> threads;
	vector<vector<pair<pattern_id_t, pattern_id_t>>> kmer_matchings(num_threads);

	vector<pair<size_t, int>> prefix_order;

	prefix_order.reserve(ht1.size());

	for (size_t i = 0; i < ht1.size(); ++i)
		prefix_order.emplace_back(ht1[i].get_size() + ht2[i].get_size(), i);

	sort(prefix_order.rbegin(), prefix_order.rend());

	atomic<int> prefix = 0;

	// Find pairs of same k-mers in both databases
	for (int i = 0; i < num_threads; ++i)
		threads.emplace_back([&, i] {
			int thread_id = i;

			vector<pair<suffix_t, pattern_id_t>> kp1, kp2;

			auto& my_res = kmer_matchings[thread_id];

			while (true)
			{
				int c_id = prefix.fetch_add(1);
				if (c_id >= ht1.size())
					break;

				int h_id = prefix_order[c_id].second;
				const auto& sht1 = ht1[h_id];
				const auto& sht2 = ht2[h_id];

				kp1.clear();
				kp2.clear();
				kp1.reserve(sht1.get_size());
				kp2.reserve(sht2.get_size());

				for (auto p = sht1.cbegin(); p != sht1.cend(); ++p)
					if (p->val != sht1.empty_value)
						kp1.emplace_back(p->key, p->val);

				for (auto p = sht2.cbegin(); p != sht2.cend(); ++p)
					if (p->val != sht2.empty_value)
						kp2.emplace_back(p->key, p->val);

				sort(kp1.begin(), kp1.end());
				sort(kp2.begin(), kp2.end());

				size_t i1 = 0;
				size_t i2 = 0;
				
				while (i1 < kp1.size() && i2 < kp2.size())
				{
					if (kp1[i1].first == kp2[i2].first)
						my_res.emplace_back(kp1[i1++].second, kp2[i2++].second);
					else if (kp1[i1].first < kp2[i2].first)
						++i1;
					else
						++i2;
				}
			}
			});

	for (auto& t : threads)
		t.join();
	threads.clear();

	// Gather kmer matchings
	vector<pair<pattern_id_t, pattern_id_t>> global_kmer_matchings;
	vector<pair<pair<pattern_id_t, pattern_id_t>, size_t>> global_kmer_matchings_counts;
	size_t n_pairs = 0;
	for (const auto& x : kmer_matchings)
		n_pairs += x.size();

	global_kmer_matchings.reserve(n_pairs);

	for (const auto& x : kmer_matchings)
		global_kmer_matchings.insert(global_kmer_matchings.end(), x.begin(), x.end());

	kmer_matchings.clear();
	kmer_matchings.shrink_to_fit();

	sort(global_kmer_matchings.begin(), global_kmer_matchings.end());

	if (!global_kmer_matchings.empty())
	{
		global_kmer_matchings_counts.emplace_back(global_kmer_matchings.front(), 1);
		for (size_t i = 1; i < global_kmer_matchings.size(); ++i)
			if (global_kmer_matchings[i] == global_kmer_matchings_counts.back().first)
				++global_kmer_matchings_counts.back().second;
			else
				global_kmer_matchings_counts.emplace_back(global_kmer_matchings[i], 1);
	}

	// Calculate matrix
	auto no_samples1 = db1.getKmersCount();
	auto no_samples2 = db2.getKmersCount();
	matrix.resize(no_samples1 + no_samples2);
//	matrix.clear();

	size_t buf_size = cacheBufferMb * 1000000 / sizeof(uint32_t);
	std::vector<uint32_t> patternsBuffer_cur(buf_size);
	std::vector<uint32_t> patternsBuffer_next(buf_size);

	vector<tuple<uint32_t, uint32_t, uint32_t*, uint32_t*>> pattern_ptrs_cur;
	vector<tuple<uint32_t, uint32_t, uint32_t*, uint32_t*>> pattern_ptrs_next;

	barrier bar_patterns(num_threads + 1);

	thread thr_pattern_decompress([&]{
		auto& patterns1 = db1.getPatterns();
		auto& patterns2 = db2.getPatterns();
		size_t buf_pos = 0;

		for (size_t i = 0; i < global_kmer_matchings_counts.size(); ++i)
		{
			auto pid1 = global_kmer_matchings_counts[i].first.first;
			auto pid2 = global_kmer_matchings_counts[i].first.second;

			const auto &ps1 = patterns1[pid1];
			const auto &ps2 = patterns2[pid2];

			if (buf_pos + ps1.get_num_samples() + ps1.get_num_samples() >= buf_size)
			{
				bar_patterns.arrive_and_wait();
				bar_patterns.arrive_and_wait();
				buf_pos = 0;
//				pattern_ptr
			}

			uint32_t* buf_ptr1 = patternsBuffer_next.data() + buf_pos;
			uint32_t* buf_ptr2 = patternsBuffer_next.data() + buf_pos + ps1.get_num_samples();
			decode_pattern_samples(patterns1, pid1, buf_ptr1);
			decode_pattern_samples(patterns1, pid1, buf_ptr2);

			pattern_ptrs_cur.emplace_back(ps1.get_num_samples(), ps2.get_num_samples(), buf_ptr1, buf_ptr2);
		}
		});

	vector<thread> thr_workers;
	thr_workers.reserve(num_threads);

	for (int i = 0; i < num_threads; ++i)
		thr_workers.emplace_back([&] {
			while (true)
			{
				bar_patterns.arrive_and_wait();

				if (patternsBuffer_cur.empty())
					break;

			}
		});


	atomic<size_t> km_idx = 0;
*/
}

// *****************************************************************************************
//
void SimilarityCalculator::db2db_sp(PrefixKmerDb& db1, PrefixKmerDb& db2, SparseMatrix<uint32_t>& matrix) 
{
//	auto& ht1 = db1.getHashtables();
//	auto& ht2 = db2.getHashtables();
	auto& sk1 = db1.getSuffixKmers();
	auto& sk2 = db2.getSuffixKmers();

	size_t samples_count1 = db1.getSamplesCount();
	size_t samples_count2 = db2.getSamplesCount();

	matrix.resize(samples_count1 + samples_count2);

	vector<thread> threads;
	vector<vector<pair<pattern_id_t, pattern_id_t>>> kmer_matchings(num_threads);

	vector<pair<size_t, int>> prefix_order;

//	prefix_order.reserve(ht1.size());
	prefix_order.reserve(sk1.size());

//	for (size_t i = 0; i < ht1.size(); ++i)
//		prefix_order.emplace_back(ht1[i].get_size() + ht2[i].get_size(), i);
	for (size_t i = 0; i < sk1.size(); ++i)
		prefix_order.emplace_back(sk1[i].size() + sk2[i].size(), i);

//	sort(prefix_order.rbegin(), prefix_order.rend());
	refresh::sort::pdqsort_branchless(prefix_order.rbegin(), prefix_order.rend());

	atomic<int> prefix = 0;

//	LOG_NORMAL << "\nK-mer comparing: " + to_string(0) + " of " + to_string(ht1.size()) << "\r";

	// Find pairs of same k-mers in both databases
	for (int i = 0; i < num_threads; ++i)
		threads.emplace_back([&, i] {
		int thread_id = i;

		vector<pair<suffix_t, pattern_id_t>> kp1, kp2;

		auto& my_res = kmer_matchings[thread_id];

		while (true)
		{
			int c_id = prefix.fetch_add(1);
//			if (c_id >= ht1.size())
			if (c_id >= (int)sk1.size())
				break;

			int h_id = prefix_order[c_id].second;
//			const auto& sht1 = ht1[h_id];
//			const auto& sht2 = ht2[h_id];
			auto& lsk1 = sk1[h_id];
			auto& lsk2 = sk2[h_id];

//			sort(lsk1.begin(), lsk1.end());
//			sort(lsk2.begin(), lsk2.end());
			refresh::sort::pdqsort_branchless(lsk1.begin(), lsk1.end());
			refresh::sort::pdqsort_branchless(lsk2.begin(), lsk2.end());

			size_t i1 = 0;
			size_t i2 = 0;

			while (i1 < lsk1.size() && i2 < lsk2.size())
			{
				if (lsk1[i1].first == lsk2[i2].first)
					my_res.emplace_back(lsk1[i1++].second, lsk2[i2++].second);
				else if (lsk1[i1].first < lsk2[i2].first)
					++i1;
				else
					++i2;
			}

//			if (sht1.get_size() < sht2.get_size())
/* {
				for (auto p = sht1.cbegin(); p != sht1.cend(); ++p)
					if (p->val != sht1.empty_value)
					{
						auto q = sht2.find(p->key);
						if (q != nullptr)
							my_res.emplace_back(p->val, *q);
					}
			}*/
/*			else
			{
				for (auto q = sht2.cbegin(); q != sht2.cend(); ++q)
					if (q->val != sht2.empty_value)
					{
						auto p = sht1.find(q->key);
						if (p != nullptr)
							my_res.emplace_back(*p, q->val);
					}
			}*/

		}
			});

	for (auto& t : threads)
		t.join();

//	LOG_NORMAL << "K-mer comparing: " + to_string(ht1.size()) + " of " + to_string(ht1.size()) << "\n";

	threads.clear();

	// Gather kmer matchings
	vector<pair<pattern_id_t, pattern_id_t>> global_kmer_matchings;
	vector<pair<pair<pattern_id_t, pattern_id_t>, size_t>> global_kmer_matchings_counts;
	size_t n_pairs = 0;
	for (const auto& x : kmer_matchings)
		n_pairs += x.size();

	LOG_DEBUG << "No. of matched pairs: " << n_pairs << endl;

	global_kmer_matchings.reserve(n_pairs);

	for (const auto& x : kmer_matchings)
		global_kmer_matchings.insert(global_kmer_matchings.end(), x.begin(), x.end());

	kmer_matchings.clear();
	kmer_matchings.shrink_to_fit();

//	sort(global_kmer_matchings.begin(), global_kmer_matchings.end());
//	ParallelSort(global_kmer_matchings.data(), global_kmer_matchings.size(), num_threads, &atp);
	refresh::sort::pdqsort_branchless_tp(num_threads, global_kmer_matchings.begin(), global_kmer_matchings.end(), atp);

	if (!global_kmer_matchings.empty())
	{
		global_kmer_matchings_counts.emplace_back(global_kmer_matchings.front(), 1);
		for (size_t i = 1; i < global_kmer_matchings.size(); ++i)
			if (global_kmer_matchings[i] == global_kmer_matchings_counts.back().first)
				++global_kmer_matchings_counts.back().second;
			else
				global_kmer_matchings_counts.emplace_back(global_kmer_matchings[i], 1);
	}

	// Calculate matrix
	auto no_samples1 = db1.getSamplesCount();
	auto no_samples2 = db2.getSamplesCount();
//	matrix.resize(no_samples1 + no_samples2);
	matrix.resize(no_samples1);
	//	matrix.clear();

	size_t buf_size = max<size_t>(2*(no_samples1 + no_samples2), cacheBufferMb * 1000000) / sizeof(uint32_t);
	std::vector<uint32_t> patternsBuffer_cur(buf_size);
	std::vector<uint32_t> patternsBuffer_next(buf_size);

	vector<tuple<uint32_t, uint32_t, uint32_t*, uint32_t*>> pattern_ptrs_cur;
	vector<tuple<uint32_t, uint32_t, uint32_t*, uint32_t*>> pattern_ptrs_next;

	auto t1 = std::chrono::high_resolution_clock::now();

	barrier bar_patterns(num_threads + 1);

	size_t percent = 0;

	thread thr_pattern_decompress([&] {
		auto& patterns1 = db1.getPatterns();
		auto& patterns2 = db2.getPatterns();
		size_t buf_pos = 0;

		LOG_DEBUG << "Pattern pairs: " << 0 << " of " << global_kmer_matchings_counts.size() << "\r" << std::flush;

		for (size_t i = 0; i < global_kmer_matchings_counts.size(); ++i)
		{
			if (i * 100 / global_kmer_matchings_counts.size() >= percent) {
				LOG_NORMAL << "\r" << percent << "%" << std::flush;
				++percent;
			}

			auto pid1 = global_kmer_matchings_counts[i].first.first;
			auto pid2 = global_kmer_matchings_counts[i].first.second;

			const auto& ps1 = patterns1[pid1];
			const auto& ps2 = patterns2[pid2];

			if (buf_pos + ps1.get_num_samples() + ps2.get_num_samples() >= buf_size)
			{
				bar_patterns.arrive_and_wait();
				swap(pattern_ptrs_cur, pattern_ptrs_next);
				swap(patternsBuffer_cur, patternsBuffer_next);
				bar_patterns.arrive_and_wait();

				pattern_ptrs_next.clear();

				buf_pos = 0;

				LOG_DEBUG << "Pattern pairs: " << i << " of " << global_kmer_matchings_counts.size() << "\r" << std::flush;	
			}

			uint32_t* buf_ptr1 = patternsBuffer_next.data() + buf_pos;
			buf_pos += ps1.get_num_samples();
			uint32_t* buf_ptr2 = patternsBuffer_next.data() + buf_pos;
			buf_pos += ps2.get_num_samples();

			decode_pattern_samples(patterns1, pid1, buf_ptr1);
			decode_pattern_samples(patterns2, pid2, buf_ptr2);

			pattern_ptrs_next.emplace_back(ps1.get_num_samples(), ps2.get_num_samples(), buf_ptr1, buf_ptr2);
		}

		bar_patterns.arrive_and_wait();
		swap(pattern_ptrs_cur, pattern_ptrs_next);
		swap(patternsBuffer_cur, patternsBuffer_next);
		bar_patterns.arrive_and_wait();

		LOG_DEBUG << "Pattern: " << global_kmer_matchings_counts.size() << 
			" of " << global_kmer_matchings_counts.size() << "\r" << std::flush;
		
		if (!patternsBuffer_cur.empty())
		{
			bar_patterns.arrive_and_wait();
			patternsBuffer_cur.clear();
			bar_patterns.arrive_and_wait();
		}
		});

	vector<thread> thr_workers;
	thr_workers.reserve(num_threads);

	for (int i = 0; i < num_threads; ++i)
		thr_workers.emplace_back([&, i] {
			uint32_t thread_id = i;
			uint32_t my_mod = thread_id % num_threads;
			uint32_t my_mod1 = thread_id % num_threads;
			uint32_t my_mod2 = 2 * num_threads - 1 - my_mod1; 
			uint32_t pid = 0;

			vector<uint32_t> my_js;
			
			const uint32_t pf_skip = 32;

			bar_patterns.arrive_and_wait();
			bar_patterns.arrive_and_wait();

			while (true)
			{
				if (patternsBuffer_cur.empty())
					break;

				for (const auto &pat : pattern_ptrs_cur)
				{
					uint32_t to_add = global_kmer_matchings_counts[pid++].second; 
					uint32_t* pat_data1 = get<2>(pat);

#if 0
					for (uint32_t j = 0; j < get<0>(pat); ++j)
					{
						uint32_t row_id = pat_data1[j];

						if (row_id % num_threads != my_mod)
							continue;

						uint32_t* pat_data2 = get<3>(pat);
						auto& row_data = matrix[row_id];

						for (uint32_t k = 0; k < get<1>(pat); ++k)
							row_data[pat_data2[k]] += to_add;
					}
#else
					my_js.clear();

					uint32_t n_data1 = get<0>(pat);
					uint32_t n_data2 = get<1>(pat);
					uint32_t* pat_data2 = get<3>(pat);

					for (uint32_t j = 0; j < n_data1; ++j)
					{
						uint32_t row_id_nt = pat_data1[j] % (2 * num_threads);
						if (row_id_nt == my_mod1 || row_id_nt == my_mod2)
							my_js.emplace_back(j);
					}

					if (n_data2 < pf_skip)
					{
						for (uint32_t jj = 0; jj < my_js.size(); ++jj)
						{
							uint32_t j = my_js[jj];
							uint32_t row_id = pat_data1[j];
							auto& row_data = matrix[row_id];

							if (jj + 1 < my_js.size())
							{
								uint32_t j_next = my_js[jj+1];
								uint32_t row_id_next = pat_data1[j_next];
								auto& row_data_next = matrix[row_id_next];
								for (uint32_t k = 0; k < n_data2; ++k)
									row_data_next.prefetch(pat_data2[k]);
							}

							for (uint32_t k = 0; k < n_data2; ++k)
								row_data[pat_data2[k]] += to_add;
						}
					}
					else
					{
						for (uint32_t jj = 0; jj < my_js.size(); ++jj)
						{
							uint32_t j = my_js[jj];
							uint32_t row_id = pat_data1[j];
							auto& row_data = matrix[row_id];

							uint32_t row_id_next = (jj + 1 < my_js.size()) ? pat_data1[my_js[jj + 1]] : 0;

							for (uint32_t k = 0; k < n_data2 - pf_skip; ++k)
							{
								row_data.prefetch(pat_data2[k + pf_skip]);
								row_data[pat_data2[k]] += to_add;
							}

							for (uint32_t k = n_data2 - pf_skip; k < n_data2; ++k)
							{
								matrix[row_id_next].prefetch(pat_data2[k - (n_data2 - pf_skip)]);
								row_data[pat_data2[k]] += to_add;
							}
						}
					}
#endif
				}

				bar_patterns.arrive_and_wait();
				bar_patterns.arrive_and_wait();
			}
			});


	thr_pattern_decompress.join();
	for (auto& t : thr_workers)
		t.join();

	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dt1 = t2 - t1;

	LOG_VERBOSE << "Db2Db Sparse computation time: " << dt1.count() << endl;
}
