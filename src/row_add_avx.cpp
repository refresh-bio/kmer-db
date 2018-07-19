/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "row_add.h"
#include <emmintrin.h>


// *****************************************************************************************
//
#ifdef NO_AVX2
void row_add(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add, bool avx2_present) {
	row_add_avx(row, src_ids, num_elems, to_add);
}
#endif


// *****************************************************************************************
//
void row_add_avx(uint32_t *row, uint32_t *src_ids, uint32_t num_elems, uint32_t to_add)
{
	__m128i _to_add = _mm_set1_epi32((int)to_add);
	auto p = src_ids;

	int j;

	if (num_elems % 32 >= 16)
	{
		j = -16;
		goto inner_start;
	}

	for (j = 0; j + 32 <= num_elems; j += 32)
	{
		if (*p + 15 == *(p + 15))
		{
			auto _q = (__m128i*) (row + *p);

			_mm_storeu_si128(_q, _mm_add_epi32(_mm_loadu_si128(_q), _to_add));
			_mm_storeu_si128(_q + 1, _mm_add_epi32(_mm_loadu_si128(_q + 1), _to_add));
			_mm_storeu_si128(_q + 2, _mm_add_epi32(_mm_loadu_si128(_q + 2), _to_add));
			_mm_storeu_si128(_q + 3, _mm_add_epi32(_mm_loadu_si128(_q + 3), _to_add));

			p += 16;
		}
		else
		{
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
		}

	inner_start:
		if (*p + 15 == *(p + 15))
		{
			auto _q = (__m128i*) (row + *p);

			_mm_storeu_si128(_q, _mm_add_epi32(_mm_loadu_si128(_q), _to_add));
			_mm_storeu_si128(_q + 1, _mm_add_epi32(_mm_loadu_si128(_q + 1), _to_add));
			_mm_storeu_si128(_q + 2, _mm_add_epi32(_mm_loadu_si128(_q + 2), _to_add));
			_mm_storeu_si128(_q + 3, _mm_add_epi32(_mm_loadu_si128(_q + 3), _to_add));

			p += 16;
		}
		else
		{
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
			row[*p++] += to_add;
		}
	}
	num_elems -= j;

	switch (num_elems % 16)
	{
	case 15:	row[*p++] += to_add;
	case 14:	row[*p++] += to_add;
	case 13:	row[*p++] += to_add;
	case 12:	row[*p++] += to_add;
	case 11:	row[*p++] += to_add;
	case 10:	row[*p++] += to_add;
	case 9:		row[*p++] += to_add;
	case 8:		row[*p++] += to_add;
	case 7:		row[*p++] += to_add;
	case 6:		row[*p++] += to_add;
	case 5:		row[*p++] += to_add;
	case 4:		row[*p++] += to_add;
	case 3:		row[*p++] += to_add;
	case 2:		row[*p++] += to_add;
	case 1:		row[*p++] += to_add;
	}
}
