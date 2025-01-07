#pragma once
#include "prefix_kmer_db.h"

class SimilarityCalculator {
public:
	SimilarityCalculator(int _num_threads, size_t cacheBufferMb);

	void all2all(PrefixKmerDb& db, LowerTriangularMatrix<uint32_t>& matrix) const;
	void all2all_sp(PrefixKmerDb& db, SparseMatrix<uint32_t>& matrix, CBubbleHelper& bubbles) const;

	template <bool parallel = true>
	void one2all(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<uint32_t>& similarities) const;
	void one2all_sp(const PrefixKmerDb& db, const kmer_t* kmers, size_t kmersCount, std::vector<std::pair<sample_id_t,num_kmers_t>>& similarities) const;

	void db2db(const PrefixKmerDb& db1, const PrefixKmerDb& db2, LowerTriangularMatrix<uint32_t>& matrix) const;
	void db2db_sp(PrefixKmerDb& db1, PrefixKmerDb& db2, SparseMatrix<uint32_t>& matrix, CBubbleHelper& bubbles) const;

protected:

	static const int PREFETCH_DIST = 48;

	int num_threads;

	size_t cacheBufferMb;

	mutable refresh::active_thread_pool atp;

	bool avx2_present;

	template<bool use_prefetch = true>
	int decode_pattern_samples_prefetch(const vector<pattern_t>& patterns, int pid, uint32_t* samples) const
	{
		const auto& pattern = patterns[pid];
		int num_samples = pattern.get_num_samples();

		uint32_t* out = samples + pattern.get_num_samples(); // start from the end

		int64_t current_id = pid;
		while (current_id >= 0) {
			const auto& cur = patterns[current_id];
			auto parent_id = cur.get_parent_id();

			if (parent_id >= 0)
#ifdef WIN32
				_mm_prefetch((const char*)(patterns.data() + parent_id), _MM_HINT_T0);
#else
				__builtin_prefetch(patterns.data() + parent_id);
#endif

			out -= cur.get_num_local_samples();
			cur.decodeSamples(out);

			current_id = parent_id;
		}

		return num_samples;
	}

	int decode_pattern_samples(const vector<pattern_t>& patterns, int pid, uint32_t* samples) const
	{
		const auto& pattern = patterns[pid];
		int num_samples = pattern.get_num_samples();

		uint32_t* out = samples + pattern.get_num_samples(); // start from the end

		int64_t current_id = pid;
		while (current_id >= 0) {
			const auto& cur = patterns[current_id];

			out -= cur.get_num_local_samples();
			cur.decodeSamples(out);

			current_id = cur.get_parent_id();
		}

		return num_samples;
	}
};