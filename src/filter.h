#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "types.h"

#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <set>


class AbstractFilter {
public:
	virtual std::unique_ptr<AbstractFilter> clone() const = 0;
	virtual void setParams(uint64_t kmerLength) { this->kmerLength = kmerLength; }
	virtual ~AbstractFilter() {}

protected:
	uint64_t kmerLength;
};

// *****************************************************************************************
//

class SetFilter : public AbstractFilter {
public:
	SetFilter(size_t initialSize) {
		//kmers.reserve_for_additional(initialSize);
	}
	void add(kmer_t kmer) { kmers.insert(kmer); }
	bool operator()(kmer_t kmer) { return kmers.find(kmer) != kmers.end(); }
	
	std::unique_ptr<AbstractFilter> clone() const override {
		return std::unique_ptr<AbstractFilter>(new SetFilter(*this));
		return nullptr;
	}

protected:
	std::set<kmer_t> kmers;
};


// *****************************************************************************************
//
class MinHashFilter : public AbstractFilter {
public:

	uint64_t getLength() const { return kmer_length; }
	double getFilterValue() const { return filterValue; }

	MinHashFilter(double filterValue, uint64_t kmer_length) : kmer_length(kmer_length), filterValue(filterValue) {
		if (filterValue <= 0.99999) {
			threshold = (uint64_t)((double)std::numeric_limits<uint64_t>::max() * filterValue);
			mode = THRESHOLD_BASED;
		}
		else {
			mode = PASS_ALL;
		}

		setParams(kmer_length);
	}

	bool operator()(kmer_t kmer) const {
		if (mode == PASS_ALL) {
			return true;
		}

		uint64_t h = hash(kmer);
		return (h < threshold);
	}
	
	void setParams(uint64_t kmer_length) override {
		AbstractFilter::setParams(kmer_length);
		this->k_div_4 = (uint64_t)ceil((double)kmer_length / 4);
		this->c42_xor_k_div_4 = 42 ^ k_div_4;
	}

	std::unique_ptr<AbstractFilter> clone() const override {
		return unique_ptr<AbstractFilter>(new MinHashFilter(*this));
	}


protected:

	enum Mode { PASS_ALL, THRESHOLD_BASED, SKETCH_BASED };

	uint64_t kmer_length; 
	double filterValue;
	
	Mode mode;
	uint64_t threshold;
	uint64_t k_div_4;
	uint64_t c42_xor_k_div_4;

	
	FORCE_INLINE uint64_t fmix64(uint64_t k) const
	{
		k ^= k >> 33;
		k *= 0xff51afd7ed558ccdull;
		k ^= k >> 33;
		k *= 0xc4ceb9fe1a85ec53ull;
		k ^= k >> 33;

		return k;
	}


	uint64_t rotl64(uint64_t x, int32_t offset) const
	{
#ifdef WIN32
		return _rotl64(x, offset);
#else
		return (x << offset) | (x >> (64 - offset));
#endif
	}

	uint64_t hash(kmer_t kmer) const {
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


};

// *****************************************************************************************
//
class SketchFilter : public MinHashFilter {
	


public:
	
	// *****************************************************************************************
	//
	SketchFilter(double filterValue, uint64_t kmer_length) : MinHashFilter(filterValue, kmer_length)
	{
		if (filterValue <= 0.99999) {
			threshold = (uint64_t)((double)std::numeric_limits<uint64_t>::max() * filterValue);
			mode = THRESHOLD_BASED;
		}
		else if (filterValue >= 2.0) {
			mode = SKETCH_BASED;
			kmer_heap.resize((int)std::round(filterValue));
		}
		
		setParams(kmer_length);
	}

	
	// *****************************************************************************************
	//
	virtual void initialize(size_t maxKmersCount) {
		if (mode == THRESHOLD_BASED || mode == PASS_ALL) {
			kmers.resize(maxKmersCount);
			kmersCount = 0;
		}
	}

	// *****************************************************************************************
	//
	virtual std::vector<kmer_t> finalize() {
		if (mode == SKETCH_BASED) {
			kmers.resize(kmer_heap.size());
			std::transform(kmer_heap.getData().begin(), kmer_heap.getData().end(), kmers.begin(),
				[](hashed_kmer_t& a)->kmer_t { return a.kmer; });
			
		}
		
		kmers.resize(kmersCount);
		return std::move(kmers);
	}

	// *****************************************************************************************
	//
	virtual void add(kmer_t kmer) const {
		// perform hashing only when needed
		if (mode == PASS_ALL) {
			kmers[kmersCount++] = kmer;
			return;
		}

		uint64_t h = hash(kmer);

		if (mode == THRESHOLD_BASED) {
			// threshold-based filtering
			if (h < threshold) {
				kmers[kmersCount++] = kmer;
			}
		}
		else {
			kmer_heap.push(hashed_kmer_t(h, kmer));
		}
	}

	
private:

	struct hashed_kmer_t {
		uint64_t hash;
		kmer_t kmer;
		
		hashed_kmer_t() {}
		hashed_kmer_t(uint64_t hash, kmer_t kmer) : hash(hash), kmer(kmer) {}
		bool operator<(hashed_kmer_t& a) const { return this->hash < a.hash; }
	};

	mutable size_t kmersCount = 0;
	mutable std::vector<kmer_t> kmers;
	mutable heap<hashed_kmer_t> kmer_heap;
};