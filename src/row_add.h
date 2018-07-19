#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <cstdint>

// uncomment this to disable AVX2 compilation
//#define NO_AVX2

void row_add(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add, bool avx2_present);

void row_add_avx(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);
void row_add_avx2(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);

