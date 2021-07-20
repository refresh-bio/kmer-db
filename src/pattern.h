#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "types.h"
#include "elias_gamma.h"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <utility>



class pattern_minimal_t {
private:
	int64_t num_kmers;			// number of kmers with this pattern
	int64_t parent_id;			// parrent pattern id

	sample_id_t num_samples;	// number of samples in current node and its parents (cannot be larger than sample id)	
/*
	pattern_minimal_t(pattern_t &v, int64_t parent_id, uint64_t num_kmers) :
		num_kmers(num_kmers), parent_id(-1), num_samples(v.num_samples + 1)
	{
		if (v.num_samples > 0) {
			v.is_parent = true;
			this->parent_id = parent_id;
		}
	}
*/
};



// *****************************************************************************************
//
class pattern_t {
private:
	int64_t num_kmers;			// number of kmers with this pattern
	int64_t parent_id;			// parrent pattern id
	sample_id_t num_samples;			// number of samples in current node and its parents (cannot be larger than sample id)	
	sample_id_t num_local_samples;		// number of samples in current node (cannot be larger than sample id)	
	sample_id_t last_sample_id;
	uint32_t num_bits;
	
	bool is_parent;			// tells whether the node is parent of some other

	uint64_t* data;				// array of samples id (Elias-Gamma encoded)

	static CEliasGamma elias;

public:
	
	// *****************************************************************************************
	//
	uint64_t* get_data() const { return data; }

	int64_t get_num_kmers() const { return num_kmers; }

	void set_num_kmers(int64_t v) { num_kmers = v; }

	void add_num_kmers(int64_t v) { num_kmers += v; }

	sample_id_t get_num_samples() const { return num_samples; }

	sample_id_t get_num_local_samples() const { return num_local_samples; }

	size_t get_num_bits() const { return num_bits; }

	bool get_is_parrent() const { return is_parent; }

	int64_t get_parent_id() const { return parent_id; }

	size_t get_data_bytes() const {
		return (num_bits == 0) ? 0 : ((num_bits + 127) / 128) * 16;
	}

	size_t get_bytes(void) const {
		return sizeof(pattern_t) + get_data_bytes();
	}

	// *****************************************************************************************
	//
	pattern_t() :
		num_kmers(0), parent_id(-1), num_samples(0), num_local_samples(0), 
		last_sample_id(0), num_bits(0), is_parent(false), data(nullptr)
	{
	}

	// *****************************************************************************************
	//
	pattern_t(sample_id_t x, uint64_t num_kmers) :
		num_kmers(num_kmers), parent_id(-1), num_samples(1), num_local_samples(1), 
		last_sample_id(x), num_bits(0), is_parent(false), data(nullptr)
		
	{
	}

	// *****************************************************************************************
	//
	pattern_t(pattern_t &v, int64_t parent_id, sample_id_t x, uint64_t num_kmers) :
		num_kmers(num_kmers), parent_id(-1), num_samples(v.num_samples + 1), 
		num_local_samples(1), last_sample_id(x), num_bits(0), is_parent(false), data(nullptr)
	{
		if (v.num_samples > 0)  {
			v.is_parent = true;
			this->parent_id = parent_id;
		}
	}


	// *****************************************************************************************
	//
/*	void from_minimal(const pattern_minimal_t& ref, sample_id_t sample_id) {
		this->num_kmers = ref.num_kmers;
		this->parent_id = ref.parent_id;
		this->num_samples = ref.num_samples;
		this->last_sample_id = sample_id;
		this->num_local_samples = 1;
	}
*/
	// *****************************************************************************************
	//
	pattern_t(const pattern_t &v) = delete;

	// *****************************************************************************************
	//
	pattern_t(pattern_t &&v)
	{
		*this = std::move(v);
	}

	// *****************************************************************************************
	//
	~pattern_t()
	{
		if (data)
		{
#ifdef USE_MALLOC
			free(data);
#else
			delete[] data;
#endif
		}
	}

	// *****************************************************************************************
	//
	pattern_t& operator=(const pattern_t &v) = delete;
	
	// *****************************************************************************************
	//
	pattern_t& operator=(pattern_t &&v) {
		is_parent = v.is_parent;
		num_samples = v.num_samples;
		num_local_samples = v.num_local_samples;
		num_kmers = v.num_kmers;
		data = v.data;
		parent_id = v.parent_id;
		num_bits = v.num_bits;
		last_sample_id = v.last_sample_id;

		v.data = nullptr;
		v.parent_id = -1;
		v.num_kmers = 0;
		v.num_local_samples = 0;
		v.num_samples = 0;
		v.num_bits = 0;
		v.last_sample_id = 0;
		
		return *this;
	}
	
	// *****************************************************************************************
	//
	void release() {
		if (data)
		{
#ifdef USE_MALLOC
			free(data);
#else
			delete[] data;
#endif
			data = nullptr;
		}
	}

	// *****************************************************************************************
	//
	void expand(const sample_id_t x)
	{
		++num_samples;
		++num_local_samples;

		uint32_t delta = x - last_sample_id;
		elias.Encode(delta, data, num_bits);
		last_sample_id = x;
	}

	// *****************************************************************************************
	//
	char* pack(char* buffer) const;

	// *****************************************************************************************
	//
	char * unpack(char* buffer);

	// *****************************************************************************************
	//
	void decodeSamples(uint32_t* out) const;
};

