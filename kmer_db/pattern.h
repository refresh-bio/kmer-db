#pragma once
#include "elias_gamma.h"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <utility>

typedef int64_t pattern_id_t;
typedef uint32_t sample_id_t;

class pattern_t {
private:

	int64_t num_kmers;			// number of kmers with this pattern
	int64_t parent_id;				// parrent pattern id

	uint64_t* data;					// tablica id próbek
	
	sample_id_t num_samples;				// liczba id próbek w wêŸle i jego rodzicach (nie mo¿e byæ wiêksza od id próbki)
	sample_id_t num_local_samples;		// liczba id próbek w wêŸle (nie mo¿e byæ wiêksza od id próbki)	
	sample_id_t last_sample_id;
	uint32_t num_bits;
	
	bool is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]

	static CEliasGamma elias;

public:

	const uint64_t* get_data() const { return data; }

	int64_t get_num_kmers() const { return num_kmers; }
	void set_num_kmers(int64_t v) { num_kmers = v; }
	void add_num_kmers(int64_t v) { num_kmers += v; }

	sample_id_t get_num_samples() const { return num_samples; }
	sample_id_t get_num_local_samples() const { return num_local_samples; }
	size_t get_num_bits() const { return num_bits; }

	bool get_is_parrent() const { return is_parent; }
	int64_t get_parent_id() { return parent_id; }
	const int64_t get_parent_id() const { return parent_id; }

	size_t get_data_bytes() const {
		return (num_bits == 0) ? 0 : ((num_bits + 127) / 128) * 16;
	}

	size_t get_bytes(void) const {
		return sizeof(pattern_t) + get_data_bytes();
	}

	pattern_t() : 
		num_samples(0), num_local_samples(0), last_sample_id(0),
		parent_id(-1), num_kmers(0), num_bits(0),
		data(nullptr), is_parent(false)
	{
	}

	pattern_t(sample_id_t x, uint64_t num_kmers) :
		num_samples(1), num_local_samples(1), last_sample_id(x),
		parent_id(-1), num_kmers(num_kmers), num_bits(0),
		data(nullptr), is_parent(false)
		
	{
	}

	pattern_t(pattern_t &v, int64_t parent_id, sample_id_t x, uint64_t num_kmers) :
		num_samples(v.num_samples + 1), num_local_samples(1), last_sample_id(x),
		parent_id(-1), num_kmers(num_kmers), num_bits(0),
		data(nullptr), is_parent(false)	
	{
		if (v.num_samples > 0)  {
			v.is_parent = true;
			this->parent_id = parent_id;
		}
	}

	pattern_t(const pattern_t &v) = delete;

	pattern_t(pattern_t &&v) 
	{
		*this = std::move(v);
	}

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

	pattern_t& operator=(const pattern_t &v) = delete;
	
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
	

	bool operator==(const pattern_t &v)
	{
		if (v.num_samples != num_samples)
			return false;
		if (v.parent_id != parent_id)
			return false;

		for (size_t i = 0; i < num_local_samples; ++i)
			if (data[i] != v.data[i])
				return false;

		return true;
	}

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

	// rozszerza listê genomów o now¹ pozycjê
	void expand(const sample_id_t x)
	{
		++num_samples;
		++num_local_samples;

		uint32_t delta = x - last_sample_id;
		elias.Encode(delta, data, num_bits);
		last_sample_id = x;
	}

	char* pack(char* buffer) const;

	char * unpack(char* buffer);

	void decodeSamples(uint32_t* out) const;
};

