#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <utility>
#include "types.h"

using namespace std;


void ParallelSort(kmer_t *arr, size_t arr_size, uint32_t n_threads);

void ParallelSort(pair<pattern_id_t, pattern_id_t*> *arr, size_t arr_size, pair<pattern_id_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads);

void ParallelSort(pair<kmer_or_pattern_t, pattern_id_t*> *arr, size_t arr_size, pair<kmer_or_pattern_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads);

void ParallelSort(pair<int, int>* arr, size_t arr_size, int n_threads);
