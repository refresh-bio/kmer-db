#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdlib>

// *****************************************************************************************
//
class IKmerFilter {
public:
	virtual bool operator()(uint64_t kmer) const = 0;
	virtual void setParams(uint64_t kmer_length) = 0;
	virtual std::unique_ptr<IKmerFilter> clone() const = 0;
};

// *****************************************************************************************
//
class MinHashFilter : public IKmerFilter {
	uint64_t rotl64(uint64 x, int32_t offset) const
	{
#ifdef WIN32
		return _rotl64(x, offset);
#else
		return (x << offset) | (x >> (64 - offset));
#endif
	}
public:
	
	// *****************************************************************************************
	//
	MinHashFilter(double fraction, size_t kmer_length) :
		threshold(fraction > 0.99999 ? std::numeric_limits<uint64_t>::max() : (uint64_t)((double)std::numeric_limits<uint64_t>::max() * fraction))
	{
		setParams(kmer_length);
	}

	// *****************************************************************************************
	//
	virtual std::unique_ptr<IKmerFilter> clone() const {
		return unique_ptr<IKmerFilter>(new MinHashFilter(*this));
	}

	// *****************************************************************************************
	//
	virtual void setParams(uint64_t kmer_length) override {
		this->k_div_4 = (uint64_t)ceil((double)kmer_length / 4);
		this->c42_xor_k_div_4 = 42 ^ k_div_4;
	}

	// *****************************************************************************************
	//
	virtual bool operator()(uint64_t kmer) const override {
		uint64_t h, h1, h2;

		h = kmer;
		h *= 0x87c37b91114253d5ull;
		h = rotl64(h, 31);
		h *= 0x4cf5ad432745937full;
		h1 = 42 ^ h; 
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
	uint64_t threshold;
	uint64_t k_div_4;
	uint64_t c42_xor_k_div_4;

	// *****************************************************************************************
	//
	FORCE_INLINE uint64_t fmix64(uint64_t k) const
	{
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccdull;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53ull;
		k ^= k >> 33;

		return k;
	}
};