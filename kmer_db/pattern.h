#pragma once
#include "elias_gamma.h"

#include <cstdint>
#include <cstddef>
#include <utility>


template<typename T>
class pattern_t {
private:

	int64_t num_kmers;			// number of kmers with this pattern
	int64_t parent_id;				// parrent pattern id

	uint64_t* data;					// tablica id próbek
	
	T num_samples;				// liczba id próbek w wêŸle i jego rodzicach (nie mo¿e byæ wiêksza od id próbki)
	T num_local_samples;		// liczba id próbek w wêŸle (nie mo¿e byæ wiêksza od id próbki)	
	T last_sample_id;
	uint32_t num_bits;
	
	bool is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]

/*	size_t round_count(size_t count) const
	{
		size_t bytes = sizeof(T) * count;

		if (bytes % 16)
			bytes = ((bytes / 16) + 1) * 16;

		return bytes / sizeof(T);
	}
*/
	static CEliasGamma elias;

public:

	const uint64_t* get_data() const { return data; }

	int64_t get_num_kmers() const { return num_kmers; }
	void set_num_kmers(int64_t v) { num_kmers = v; }

	T get_num_samples() const { return num_samples; }
	T get_num_local_samples() const { return num_local_samples; }
	size_t get_num_bits() const { return num_bits; }

	const T& operator[](size_t i) const { return data[i]; }
	T& operator[](size_t i) { return data[i]; }

	bool get_is_parrent() const { return is_parent; }
	int64_t get_parent_id() { return parent_id; }
	const int64_t get_parent_id() const { return parent_id; }

	pattern_t() : 
		num_samples(0), num_local_samples(0), last_sample_id(0),
		parent_id(-1), num_kmers(0), num_bits(0),
		data(nullptr), is_parent(false)
	{
	}

	pattern_t(T x, uint64_t num_kmers) : 
		num_samples(1), num_local_samples(1), last_sample_id(x),
		parent_id(-1), num_kmers(num_kmers), num_bits(0),
		data(nullptr), is_parent(false)
		
	{
	}

	pattern_t(pattern_t &v, int64_t parent_id, T x, uint64_t num_kmers) : 
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
			delete[] data;
		}
	}

	pattern_t<T>& operator=(const pattern_t<T> &v) = delete;
	
	pattern_t<T>& operator=(pattern_t<T> &&v) {
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
	

	bool operator==(const pattern_t<T> &v)
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
			delete[] data;
			data = nullptr;
		}
	}

	// rozszerza listê genomów o now¹ pozycjê
	void expand(const T x)
	{
		++num_samples;
		++num_local_samples;

		uint32_t delta = x - last_sample_id;
		elias.Encode(delta, data, num_bits);
		last_sample_id = x;
	}

	size_t get_data_bytes() const { 
		return (num_bits == 0) ? 0 : (num_bits + 127) / 8; 
	}

	size_t get_bytes(void) const { 
		return sizeof(pattern_t) + get_data_bytes();
	}

	
	char* pack(char* buffer) const {
		// store members
		*reinterpret_cast<decltype(num_kmers)*>(buffer) = num_kmers;
		buffer += sizeof(decltype(num_kmers));

		*reinterpret_cast<decltype(parent_id)*>(buffer) = parent_id;
		buffer += sizeof(decltype(parent_id));

		*reinterpret_cast<decltype(is_parent)*>(buffer) = is_parent;
		buffer += sizeof(decltype(is_parent));

		*reinterpret_cast<decltype(num_samples)*>(buffer) = num_samples;
		buffer += sizeof(decltype(num_samples));

		*reinterpret_cast<decltype(num_local_samples)*>(buffer) = num_local_samples;
		buffer += sizeof(decltype(num_local_samples));

		*reinterpret_cast<decltype(num_bits)*>(buffer) = num_bits;
		buffer += sizeof(decltype(num_bits));

		*reinterpret_cast<decltype(last_sample_id)*>(buffer) = last_sample_id;
		buffer += sizeof(decltype(last_sample_id));
		
		size_t data_bytes = get_data_bytes();
		if (data_bytes) {
			memcpy(buffer, reinterpret_cast<char*>(data), data_bytes);
			buffer += data_bytes;
		}

		return buffer;
	}

	char * unpack(char* buffer) {
		if (num_local_samples) {
			delete[] data;
		}
		
		num_kmers = *reinterpret_cast<decltype(num_kmers)*>(buffer);
		buffer += sizeof(decltype(num_kmers));

		parent_id = *reinterpret_cast<decltype(parent_id)*>(buffer);
		buffer += sizeof(decltype(parent_id));

		is_parent = *reinterpret_cast<decltype(is_parent)*>(buffer);
		buffer += sizeof(decltype(is_parent));

		num_samples = *reinterpret_cast<decltype(num_samples)*>(buffer);
		buffer += sizeof(decltype(num_samples));

		num_local_samples = *reinterpret_cast<decltype(num_local_samples)*>(buffer);
		buffer += sizeof(decltype(num_local_samples));
		
		num_bits = *reinterpret_cast<decltype(num_bits)*>(buffer);
		buffer += sizeof(decltype(num_bits));

		last_sample_id = *reinterpret_cast<decltype(last_sample_id)*>(buffer);
		buffer += sizeof(decltype(last_sample_id));


		size_t num_bytes = get_data_bytes();

		if (num_bytes) {
			data = new uint64_t[num_bytes / sizeof(uint64_t)]; 
			std::memcpy(reinterpret_cast<char*>(data), buffer, num_bytes);
			buffer += num_bytes;
		}

		return buffer;
	}

	void decodeSamples(uint32_t* out) const {
		if (num_local_samples) {
			out[num_local_samples - 1] = last_sample_id;
			if (num_local_samples > 1) {
				elias.Decode(data, num_bits, out);
				for (int i = num_local_samples - 2; i >= 0; --i) {
					out[i] = out[i + 1] - out[i];
				}
			}
		}
	}

};

template <typename T>
CEliasGamma pattern_t<T>::elias;