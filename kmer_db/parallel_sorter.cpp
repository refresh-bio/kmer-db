#include "parallel_sorter.h"

#ifdef USE_RADULS
#include "raduls/raduls.h"
#else
#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif
#endif


void ParallelSort(pair<pattern_id_t, pattern_id_t*> *arr, size_t arr_size, pair<pattern_id_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads)
{
#ifdef USE_RADULS
	raduls::RadixSortMSD(reinterpret_cast<uint8_t*>(arr), reinterpret_cast<uint8_t*>(tmp), arr_size, rec_size, key_size, n_threads);
#else
	auto pid_comparer = [](const std::pair<pattern_id_t, pattern_id_t*>& a, const std::pair<pattern_id_t, pattern_id_t*>& b)->bool {
		return a.first < b.first;
	};
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size, pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#else
	__gnu_parallel::sort(arr, arr + size, pid_comparer);
#endif
#endif
}
