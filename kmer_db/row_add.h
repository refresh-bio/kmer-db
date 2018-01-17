#pragma once

#include <cstdint>

void row_add_avx(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);
void row_add_avx2(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add);

