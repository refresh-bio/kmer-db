#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include <cstdint>

void row_add_avx(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);
void row_add_avx2(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);

