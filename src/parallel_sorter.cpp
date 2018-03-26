/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include "parallel_sorter.h"

#ifdef WIN32
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif

// *****************************************************************************************
//
void ParallelSort(pair<pattern_id_t, pattern_id_t*> *arr, size_t arr_size, pair<pattern_id_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads)
{
	auto pid_comparer = [](const std::pair<pattern_id_t, pattern_id_t*>& a, const std::pair<pattern_id_t, pattern_id_t*>& b)->bool {
		return a.first < b.first;
	};
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size, pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#else
	__gnu_parallel::sort(arr, arr + arr_size, pid_comparer);
#endif
}
