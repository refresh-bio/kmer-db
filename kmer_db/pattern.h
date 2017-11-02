#pragma once


template<typename T>
class pattern_t {
private:
	bool is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]

	T num_samples;				// liczba id próbek w wêŸle i jego rodzicach (nie mo¿e byæ wiêksza od id próbki)
	T num_local_samples;		// liczba id próbek w wêŸle (nie mo¿e byæ wiêksza od id próbki)	
	T* data;					// tablica id próbek

	int64_t parent_id;				// parrent pattern id

							
	size_t round_count(size_t count)
	{
		size_t bytes = sizeof(T) * count;

		if (bytes % 16)
			bytes = ((bytes / 16) + 1) * 16;

		return bytes / sizeof(T);
	}


public:

	uint64_t num_kmers;			// number of kmers with this pattern


	const T* get_data() const { return data; }

	size_t get_num_samples() const { return num_samples; }
	size_t get_num_local_samples() const { return num_local_samples; }

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
		if (v.parent != parent)
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

	size_t get_bytes(void)
	{
		size_t r;

		r = sizeof(pattern_t);
		r += sizeof(T) * round_count(num_local_samples);

		return r;
	}

	std::string toString() const {
		std::string out;
		if (parent) {
			out += parent->toString() + " | ";
		}
		
		for (int i = 0; i< num_local_samples; ++i) {
			out += std::to_string(data[i]) + ", ";
		}
		return out;
	}
};


/*
template <class T>
class pattern_t
{
public:
	uint32_t num_occ;
	subpattern_t<T> *last_subpattern;

	std::string toString() const {
		return last_subpattern->toString();
	}

	std::vector<T> toSamplesVector() const {
		std::vector<T> samples;

		auto subpattern = last_subpattern;
		samples.resize(subpattern->get_num_samples());
		int j = samples.size() - 1;

		// collect in reversed order
		do {
			for (int i = subpattern->get_num_local_samples() - 1; i >= 0; --i, --j) {
				samples[j] = (*subpattern)[i];
			}
			subpattern = subpattern->get_parent();

		} while (subpattern != nullptr);

		return samples;
	}

};
*/