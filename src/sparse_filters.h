#pragma once
#include <vector>
#include <limits>
#include <map>
#include <string>

#include "types.h"

using metric_fun_t = std::function<double(num_kmers_t, num_kmers_t, num_kmers_t, int)>;

// Filter applied on metric
struct MetricFilter {

	double bounds[2] { std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max() };

	metric_fun_t metric;

	// calculate measure and check if passed
	bool operator()(num_kmers_t common, num_kmers_t cnt1, num_kmers_t cnt2, int kmerLength) const {
		double value = metric(common, cnt1, cnt2, kmerLength);
		return (value >= bounds[0] && value <= bounds[1]);
	}
};

// Filter applied on kmer count
struct KmerFilter {
	num_kmers_t bounds[2] { 0,  std::numeric_limits<num_kmers_t>::max() };

	bool operator()(num_kmers_t n) const { return n >= bounds[0] && n <= bounds[1]; }
};


template <class T>
struct CombinedFilter {

	const std::map<std::string,MetricFilter>& metricFilters;
	const KmerFilter& kmerFilter;
	const std::vector<T>& rowCounts;
	const std::vector<T>& colCounts;
	int kmerLength;

	CombinedFilter(const std::map<std::string, MetricFilter>& metricFilters, const KmerFilter& kmerFilter, const std::vector<T>& rowCounts, const std::vector<T>& colCounts, int kmerLength) :
		metricFilters(metricFilters),
		kmerFilter(kmerFilter),
		rowCounts(rowCounts),
		colCounts(colCounts),
		kmerLength(kmerLength) {}

	bool operator()(T common, int row_id, int col_id) const {
		for (const auto& filter : metricFilters) {
			if (!filter.second(common, rowCounts[row_id], colCounts[col_id], kmerLength)) {
				return false;
			}
		}

		if (!kmerFilter(common)) {
			return false;
		}
		return true;
	}

};