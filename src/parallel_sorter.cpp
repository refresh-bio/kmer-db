/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#define USE_PDQSORT

#include "parallel_sorter.h"
#include "types.h"

#ifdef USE_PDQSORT
#include "../libs/refresh/pdqsort_par.h"
#else
#ifdef WIN32
#include <ppl.h>
#elif defined __APPLE__
#include <algorithm>
#else
#include <parallel/algorithm>
#endif
#endif

// *****************************************************************************************
//
void ParallelSort(kmer_t *arr, size_t arr_size, uint32_t max_n_threads)
{
#ifdef USE_PDQSORT
	refresh::sort::pdqsort_branchless(refresh::sort::pdqsort_adjust_threads(arr_size, max_n_threads), arr, arr + arr_size);
#else
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#elif defined __APPLE__
	std:: stable_sort(arr, arr + arr_size);
#else
	__gnu_parallel::sort(arr, arr + arr_size);
#endif
#endif
}

// *****************************************************************************************
//
void ParallelSort(pair<pattern_id_t, pattern_id_t*> *arr, size_t arr_size, pair<pattern_id_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads)
{
	auto pid_comparer = [](const std::pair<pattern_id_t, pattern_id_t*>& a, const std::pair<pattern_id_t, pattern_id_t*>& b)->bool {
		return a.first < b.first;
	};

#ifdef USE_PDQSORT
	refresh::sort::pdqsort_branchless(refresh::sort::pdqsort_adjust_threads(arr_size, n_threads), arr, arr + arr_size, pid_comparer);
#else
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size, pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#elif defined __APPLE__
	std:: stable_sort(arr, arr + arr_size, pid_comparer);
#else
	__gnu_parallel::sort(arr, arr + arr_size, pid_comparer);
#endif
#endif
}

// *****************************************************************************************
//
void ParallelSort(pair<int, int>* arr, size_t arr_size, int n_threads)
{
	auto pid_comparer = [](const std::pair<int, int>& a, const std::pair<int, int>& b)->bool {
		return a < b;
		};

#ifdef USE_PDQSORT
	refresh::sort::pdqsort_branchless(refresh::sort::pdqsort_adjust_threads(arr_size, n_threads), arr, arr + arr_size, pid_comparer);
#else
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size, pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#elif defined __APPLE__
	std::stable_sort(arr, arr + arr_size, pid_comparer);
#else
	__gnu_parallel::sort(arr, arr + arr_size, pid_comparer);
#endif
#endif
}

// *****************************************************************************************
//
void ParallelSort(pair<kmer_or_pattern_t, pattern_id_t*> *arr, size_t arr_size, pair<kmer_or_pattern_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads)
{
	auto pid_comparer = [](const std::pair<kmer_or_pattern_t, pattern_id_t*>& a, const std::pair<kmer_or_pattern_t, pattern_id_t*>& b)->bool {
		return a.first.pattern_id < b.first.pattern_id;
	};

#ifdef USE_PDQSORT
	refresh::sort::pdqsort_branchless(refresh::sort::pdqsort_adjust_threads(arr_size, n_threads), arr, arr + arr_size, pid_comparer);
#else
#ifdef WIN32
	concurrency::parallel_sort(arr, arr + arr_size, pid_comparer);
	//std::stable_sort(samplePatterns.begin(), samplePatterns.end(), pid_comparer);
#elif defined __APPLE__
	std:: stable_sort(arr, arr + arr_size, pid_comparer);
#else
	__gnu_parallel::sort(arr, arr + arr_size, pid_comparer);
#endif
#endif
}
