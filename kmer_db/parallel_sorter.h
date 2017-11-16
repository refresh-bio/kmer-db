#pragma once

#include <utility>
#include "pattern.h"

#define USE_RADULS

using namespace std;

void ParallelSort(pair<pattern_id_t, pattern_id_t*> *arr, size_t arr_size, pair<pattern_id_t, pattern_id_t*> *tmp, int rec_size, int key_size, int n_threads);

