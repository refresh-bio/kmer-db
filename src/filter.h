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

	IKmerFilter(double filterValue, uint64_t kmer_length) : kmer_length(kmer_length), filterValue(filterValue) {}

	virtual void setParams(uint64_t kmer_length) = 0;
	uint64_t getLength() const { return kmer_length; }
	double getFilterValue() const { return filterValue; }
	virtual bool operator()(uint64_t kmer) const = 0;
	virtual std::unique_ptr<IKmerFilter> clone() const = 0;

protected:
	uint64_t kmer_length; 
	double filterValue;
};

// *****************************************************************************************
//
class MinHashFilter : public IKmerFilter {
	
	enum Mode {PASS_ALL, THRESHOLD_BASED, SKETCH_BASED};
	
	uint64_t rotl64(uint64 x, int32_t offset) const
	{
#ifdef WIN32
		return _rotl64(x, offset);
#else
		return (x << offset) | (x >> (64 - offset));
#endif
	}

	uint64_t hash(uint64_t kmer) const {
		uint64_t h, h1, h2;

		// calculate hash
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

		return h1 ^ h2; // xor as final hash
	}


public:
	
	// *****************************************************************************************
	//
	MinHashFilter(double filterValue, uint64_t kmer_length) : IKmerFilter(filterValue, kmer_length), sketchSize(0), threshold(std::numeric_limits<uint64_t>::max()), mode(PASS_ALL)
	{
		if (filterValue <= 0.99999) {
			threshold = (uint64_t)((double)std::numeric_limits<uint64_t>::max() * filterValue);
			mode = THRESHOLD_BASED;
		}
		else if (filterValue >= 2.0) {
			sketchSize = (int)std::round(filterValue);
			mode = SKETCH_BASED;
		}
		
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
		this->kmer_length = kmer_length;
		this->k_div_4 = (uint64_t)ceil((double)kmer_length / 4);
		this->c42_xor_k_div_4 = 42 ^ k_div_4;
	}

	// *****************************************************************************************
	//
	virtual bool operator()(uint64_t kmer) const override {
		// perform hashing only when needed
		if (mode == PASS_ALL) {
		//	collection.push_back(kmer);
			return true;
		}

		uint64_t h = hash(kmer);

		if (mode == THRESHOLD_BASED) {
			// threshold-based filtering
			if (h < threshold) {
			//	collection.push_back(kmer);
				return true;
			}
		}
		else {
			// TODO
			
			return false;
		}

		return false;
	}

private:
	Mode mode;
	int sketchSize;
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