#include "similarity_calculator.h"

#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <emmintrin.h>
#include <immintrin.h>
#include "instrset.h"


// *****************************************************************************************
//
SimilarityCalculator::SimilarityCalculator(int _num_threads, size_t cacheBufferMb) : 
	num_threads(_num_threads > 0 ? _num_threads : std::thread::hardware_concurrency()),
	cacheBufferMb(cacheBufferMb) {
	avx2_present = instrset_detect() >= 8;
}

// *****************************************************************************************
//
void SimilarityCalculator::operator()(PrefixKmerDb& db, LowerTriangularMatrix<uint32_t>& matrix) const 
{
	// get stuff from database
	auto& patterns = db.getPatterns();
	int samples_count = db.getSamplesCount();

	matrix.resize(samples_count);
	matrix.clear();

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


	// Decompress patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_decomp[tid] = std::thread([this, samples_count, &patterns, &sample2pattern, &v_range_patterns, &tasks_decomp_queue, &semaphore_decomp, &rawPatterns, &first_pid, &patternsBuffer, &hist_sample_ids]
		{
			while (!tasks_decomp_queue.IsCompleted())
			{
				pair<int, uint32_t> decomp_task;

				if (tasks_decomp_queue.Pop(decomp_task))
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
								_mm_prefetch((const char*)(patterns.data() + parent_id), _MM_HINT_T0);

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
			}
		});
	}

	// Sample2patterns
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_sample2patterns[tid] = std::thread([
			this, &patterns, &sample2pattern, &v_range_patterns, &tasks_sample2patterns_queue, 
				&semaphore_sample2patterns, &rawPatterns, &first_pid, &patternsBuffer, &hist_sample_ids]
		{
			while (!tasks_sample2patterns_queue.IsCompleted())
			{
				int part_id;

				if (tasks_sample2patterns_queue.Pop(part_id))
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
			}
		});
	}

	// increment array elements in threads
	for (int tid = 0; tid < num_threads; ++tid) {
		workers_matrix[tid] = std::thread([this, &patterns, &sample2pattern, &rawPatterns, &workerRanges, &matrix, &tasks_matrix_queue, &first_pid, &semaphore_matrix, &numAdditions] {
			uint64_t localAdditions = 0;

			int range_id;

			while (!tasks_matrix_queue.IsCompleted())
			{
				if (tasks_matrix_queue.Pop(range_id))
				{
					for (int id = workerRanges[range_id]; id < workerRanges[range_id + 1]; ) {
						int Si = get<0>(sample2pattern[id]);
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
			}
			numAdditions.fetch_add(localAdditions);
		});
	}

	// calculate histogram
	for (int tid = 0; tid < num_threads; ++tid)
		workers_histogram[tid] = std::thread([&hist_boundary_values, &tasks_histogram_queue, samples_count, this, &v_tmp,
			&hist_sample_ids, &semaphore_hist_first, &semaphore_hist_second] {
		int task_id;

		while (!tasks_histogram_queue.IsCompleted())
		{
			if (tasks_histogram_queue.Pop(task_id))
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
				semaphore_hist_first.waitForZero();

				int to_add = 0;
				for (int i = 0; i < task_id; ++i)
					to_add += hist_boundary_values[i];

				for (int j = 0; j < max_j; ++j)
					for (int i = first_sample; i < last_sample; ++i)
						hist_sample_ids[j][i] += to_add;
				semaphore_hist_second.dec();
			}
		}
	});

	// process all patterns in blocks determined by buffer size
	std::cout << std::endl;

	double decomp_time = 0;
	double hist_time = 0;
	double patterns_time = 0;

	int pid_to_cout = 0;
	for (int pid = 0; pid < patterns.size(); ) {
		if (pid >= pid_to_cout)
		{
			std::cout << pid << " of " << patterns.size() << "\r";
			fflush(stdout);
			pid_to_cout += 50000;
		}

		first_pid = pid;
		size_t samplesCount = 0;

		auto t1 = std::chrono::high_resolution_clock::now();

		int part_size = bufsize / no_hist_parts + 1;
		int next_boundary = 0;
		int part_id = 0;

		v_tmp.clear();
		v_tmp_int.clear();

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
		semaphore_decomp.waitForZero();

		auto t2 = std::chrono::high_resolution_clock::now();

		semaphore_hist_first.inc(num_threads);
		semaphore_hist_second.inc(num_threads);
		tasks_histogram_queue.PushRange(v_hist_ids);
		semaphore_hist_second.waitForZero();

		auto t3 = std::chrono::high_resolution_clock::now();

		// generate sample to pattern mapping
		int num_items = std::accumulate(hist_boundary_values.begin(), hist_boundary_values.end(), 0);
		sample2pattern.resize(num_items);		// resize to number of items in the current fragment

		semaphore_sample2patterns.inc(v_tmp_int.size());
		tasks_sample2patterns_queue.PushRange(v_tmp_int);
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

		for (int rid = 0; rid < no_ranges - 1; ++rid) {
			// make sure no single row of matrix is updated by multiple workers
			auto it = std::upper_bound(
				sample2pattern.begin() + currentIndex,
				sample2pattern.end(),
				*(sample2pattern.begin() + currentIndex - 1),
				[](auto x, auto y) {return get<0>(x) < get<0>(y); });	// Necessary, because we are looking for end of sample data

			size_t range = it - sample2pattern.begin();
			workerRanges[rid + 1] = range;
			currentIndex = range + workerBlock;

			if (currentIndex >= sample2pattern.size())
				break;
		}

		// this should never happen
		if (workerRanges[no_ranges] != sample2pattern.size()) {
			throw std::runtime_error("ERROR in FastKmerDb::calculateSimilarity(): Invalid ranges");
		}

		semaphore_matrix.inc(no_ranges);
		tasks_matrix_queue.PushRange(v_range_ids);
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

	cout << "Decomp   time: " << decomp_time << endl;
	cout << "Hist     time: " << hist_time << endl;
	cout << "Sample2p time: " << patterns_time << endl;

#ifdef ALL_STATS
	cout << "Number of additions:" << numAdditions << endl;
#endif
}

// *****************************************************************************************
//
void  SimilarityCalculator::operator()(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<uint32_t>& similarities) const {
	
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
	for (int i = 0; i < kmersCount; ++i) {

		if (i + prefetch_dist < kmersCount) {
			kmer_t prefetch_kmer = kmers[i + prefetch_dist];

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
	LOG_VERBOSE << "Pattern listing time: " << dt.count() << endl;

	start = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> workers(num_threads);
	for (int tid = 0; tid < num_threads; ++tid) {
		workers[tid] = std::thread([this, &patterns, samples_count, tid, &patterns2countVector, &localSimilarities, &range_boundaries]() {
			std::vector<uint32_t> samples(samples_count);

			size_t lo = range_boundaries[tid];
			size_t hi = range_boundaries[tid + 1];
			auto &my_localSimilarities = localSimilarities[tid];

			for (int id = lo; id < hi; ++id) {
				if (id + 1 < hi)
					_mm_prefetch((const char*)(patterns.data() + patterns2countVector[id + 1].first), _MM_HINT_T0);

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
	LOG_VERBOSE << "Pattern unpacking time: " << dt.count() << endl;
}