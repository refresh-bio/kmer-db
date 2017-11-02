#pragma once


template<typename T>
class subpattern_t {

	
	uint16_t is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]
	uint16_t is_child;			// informacja czy wêze³ jest zapisany jako rozszerzenie innego [true] czy te¿ wprost (lista id genomów) [false]
	T num_samples;				// liczba id próbek w wêŸle i jego rodzicach (nie mo¿e byæ wiêksza od id próbki)
	T num_local_samples;		// liczba id próbek w wêŸle (nie mo¿e byæ wiêksza od id próbki)
	T* data;					// tablica id próbek
	subpattern_t<T>* parent;		// wskaŸnik do wêz³a rodzica
	
							
	size_t round_count(size_t count)
	{
		size_t bytes = sizeof(T) * count;

		if (bytes % 16)
			bytes = ((bytes / 16) + 1) * 16;

		return bytes / sizeof(T);
	}


public:

	const T* get_data() const { return data; }

	size_t get_num_samples() const { return num_samples; }
	size_t get_num_local_samples() const { return num_local_samples; }

	const T& operator[](size_t i) const { return data[i]; }
	T& operator[](size_t i) { return data[i]; }

	bool get_is_parrent() const { return is_parent; }
	bool get_is_child() const { return is_child; }
	subpattern_t* get_parent() { return parent; }
	const subpattern_t* get_parent() const { return parent; }

	subpattern_t() : num_samples(0), num_local_samples(0), data(nullptr), is_parent(false), is_child(false), parent(nullptr)
	{
	}

	subpattern_t(T x)
	{
		is_parent = false;
		is_child = false;
		num_samples = 1;
		num_local_samples = 1;
		data = nullptr;
		parent = nullptr;

		data = new T[round_count(1)];
		data[0] = x;
	}

	subpattern_t(subpattern_t &v, T x)
	{
		num_samples = v.num_samples + 1;
		num_local_samples = 1;
		data = nullptr;

		if (v.num_samples == 0)
		{
			is_parent = false;
			is_child = false;
			v.is_parent = false;
			parent = nullptr;
		}
		else
		{
			is_parent = false;
			is_child = true;
			v.is_parent = true;
			parent = &v;
		}

		data = new T[round_count(1)];
		data[0] = x;
	}

	subpattern_t(const subpattern_t &v) = delete;
	/*{
		is_parrent = v.is_parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		no_locally_allocated = v.no_locally_allocated;
		data = nullptr;
		parrent = v.parrent;

		if (no_ids)
		{
			size_t to_copy = no_ids;

			if (parrent) {
				to_copy -= parrent->no_ids;
			}

			data = new T[no_locally_allocated]
			
			for (int i = 0; i < to_copy; ++i) {
				data[i] = v.data[i];
			}

			mem_pattern_desc += sizeof(T);
			mem_os_pattern_desc += sizeof(T) * round_count(to_alloc);
		}
		else
			data = nullptr;
	}*/

	subpattern_t(subpattern_t &&v) = delete;
	/*{
		is_parrent = v.is_parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		no_locally_allocated = v.no_locally_allocated;
		data = v.data;
		parrent = v.parrent;

		v.data = nullptr;
		v.parrent = nullptr;
		v.no_ids = 0;
		v.no_locally_allocated = 0;
	}*/ 

	~subpattern_t()
	{
		if (num_samples)
		{
			delete[] data;
		}
	}

	subpattern_t<T>& operator=(const subpattern_t<T> &v) = delete;
	subpattern_t<T>& operator=(subpattern_t<T> &&v) = delete;
	

	bool operator==(const subpattern_t<T> &v)
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

	size_t getMem(void)
	{
		size_t r;

		r = sizeof(subpattern_t);
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