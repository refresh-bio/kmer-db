#pragma once
#include <cstdint>
#include <cstddef>
#include <utility>

template<typename T>
class pattern_t {
private:

	int64_t num_kmers;			// number of kmers with this pattern
	int64_t parent_id;				// parrent pattern id

	T* data;					// tablica id próbek
	
	T num_samples;				// liczba id próbek w wêŸle i jego rodzicach (nie mo¿e byæ wiêksza od id próbki)
	T num_local_samples;		// liczba id próbek w wêŸle (nie mo¿e byæ wiêksza od id próbki)	
	
	bool is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]

	size_t round_count(size_t count) const
	{
		size_t bytes = sizeof(T) * count;

		if (bytes % 16)
			bytes = ((bytes / 16) + 1) * 16;

		return bytes / sizeof(T);
	}


public:

	const T* get_data() const { return data; }

	int64_t get_num_kmers() const { return num_kmers; }
	void set_num_kmers(int64_t v) { num_kmers = v; }

	T get_num_samples() const { return num_samples; }
	T get_num_local_samples() const { return num_local_samples; }

	const T& operator[](size_t i) const { return data[i]; }
	T& operator[](size_t i) { return data[i]; }

	bool get_is_parrent() const { return is_parent; }
	int64_t get_parent_id() { return parent_id; }
	const int64_t get_parent_id() const { return parent_id; }

	pattern_t() : num_samples(0), num_local_samples(0), data(nullptr), is_parent(false), parent_id(-1), num_kmers(0)
	{
	}

	pattern_t(T x, uint64_t num_kmers)
	{
		is_parent = false;
		num_samples = 1;
		num_local_samples = 1;
		data = nullptr;
		parent_id = -1;
		this->num_kmers = num_kmers;

		data = new T[round_count(1)];
		data[0] = x;
	}

	pattern_t(pattern_t &v, int64_t parent_id, T x, uint64_t num_kmers)
	{
		num_samples = v.num_samples + 1;
		num_local_samples = 1;
		data = nullptr;
		this->num_kmers = num_kmers;

		if (v.num_samples == 0)
		{
			is_parent = false;
			v.is_parent = false;
			this->parent_id = -1;
		}
		else
		{
			is_parent = false;
			v.is_parent = true;
			this->parent_id = parent_id;
		}

		data = new T[round_count(1)];
		data[0] = x;
	}

	pattern_t(const pattern_t &v) = delete;

	pattern_t(pattern_t &&v) 
	{
		*this = std::move(v);
	}

	~pattern_t()
	{
		if (num_samples)
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

		v.data = nullptr;
		v.parent_id = -1;
		v.num_kmers = 0;
		v.num_local_samples = 0;
		v.num_samples = 0;

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

	// rozszerza listê genomów o now¹ pozycjê
	void expand(const T x)
	{
		++num_samples;
		
		// if no space left - reallocate and copy existing ids 
		if (num_local_samples + 1 > round_count(num_local_samples)) {
			T* old_data = data;
		
			data = new T[round_count(num_local_samples + 1)];
			
			for (int i = 0; i < num_local_samples; ++i) {
				data[i] = old_data[i];
			}
			
			delete[] old_data;
		}

		++num_local_samples;
		data[num_local_samples - 1] = x;
	}

	size_t get_bytes(void) const
	{
		size_t r;

		r = sizeof(pattern_t);
		r += sizeof(T) * round_count(num_local_samples);

		return r;
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
		
		T* samples = reinterpret_cast<T*>(buffer);
		for (int i = 0; i < num_local_samples; ++i) {
			samples[i] = data[i];
		}

		buffer += num_local_samples * sizeof(T);

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
		
		data = new T[round_count(num_local_samples)];

		T* samples = reinterpret_cast<T*>(buffer);
		for (int i = 0; i < num_local_samples; ++i) {
			data[i] = samples[i];
		}

		buffer += num_local_samples * sizeof(T);

		return buffer;
	}

};

