#pragma once

#include <cstdint>

#define KMER_MSB (1ULL << 63)

#define SUFFIX_BITS 32

#define SUFFIX_MASK ((1ULL << SUFFIX_BITS) - 1)
#define PREFIX_MASK (~SUFFIX_MASK)

typedef uint64_t kmer_t;
typedef uint32_t suffix_t;

typedef int32_t pattern_id_t;
typedef uint32_t sample_id_t;

typedef uint32_t num_kmers_t;

union kmer_or_pattern_t {
	kmer_t kmer;
	pattern_id_t pattern_id;
};

#define GET_PREFIX(kmer)			((kmer) & PREFIX_MASK)
#define GET_SUFFIX(kmer)			(static_cast<suffix_t>((kmer) & SUFFIX_MASK))
#define GET_PREFIX_SHIFTED(kmer)	((kmer) >> SUFFIX_BITS)


