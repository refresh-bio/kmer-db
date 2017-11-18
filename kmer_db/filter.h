#pragma once

#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdlib>

class IKmerFilter {
public:
	virtual bool operator()(uint64_t kmer) = 0;
};


class NullFilter : public IKmerFilter {
public:
	virtual bool operator()(uint64_t kmer) override {
		return true;
	}
};

class MinHashFilter : public IKmerFilter {
public:
	MinHashFilter(double fraction, size_t kmer_length) :
		threshold(std::numeric_limits<uint64_t>::max() * std::min(fraction, 1.0)),
		k_div_4(std::ceil(kmer_length / 4)),
		c42_xor_k_div_4(42 ^ k_div_4)
	{}

	virtual bool operator()(uint64_t kmer) override {
		uint64_t h, h1, h2;

		h = kmer;
		h *= 0x87c37b91114253d5ull;
		h = _rotl64(h, 31);
		h *= 0x4cf5ad432745937f;
		h1 ^= 42 ^ h;
		h1 ^= k_div_4; //ceil(k / 4);
		h2 = c42_xor_k_div_4; // 42 ^ ceil(k / 4);
		h1 += h2;
		h2 += h1;
		h1 = fmix64(h1);
		h2 = fmix64(h2);
		h1 += h2;
		h2 += h1;

		// make xor
		return (h1 ^ h2) < threshold;
	}

private:
	const uint64_t threshold;
	const uint64_t k_div_4;
	const uint64_t c42_xor_k_div_4;

	FORCE_INLINE uint64_t fmix64(uint64_t k)
	{
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccd;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53;
		k ^= k >> 33;

		return k;
	}
};