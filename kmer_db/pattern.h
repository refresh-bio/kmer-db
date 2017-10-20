#pragma once


template<typename T>
class subpattern_t {

	bool is_parent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]
	bool is_child;				// informacja czy wêze³ jest zapisany jako rozszerzenie innego [true] czy te¿ wprost (lista id genomów) [false]
	size_t num_samples;				// liczba id genomów w wêŸle
	size_t num_locally_allocated;	// liczba slotów zaalokowanych w wêŸle
	T* data;					// tablica id genomów
	subpattern_t<T>* parent;		// wskaŸnik do wêz³a rodzica

								
	size_t round_count(size_t count)
	{
		size_t bytes = sizeof(T) * count;

		if (bytes % 16)
			bytes = ((bytes / 16) + 1) * 16;

		return bytes / sizeof(T);
	}


public:

	size_t num_local_samples() const {
		return num_samples - (parent != nullptr ? parent->num_samples : 0);
	}

	const T& operator[](size_t i) const { return data[i]; }
	T& operator[](size_t i) { return data[i]; }

	bool get_is_parrent() const { return is_parent; }
	bool get_is_child() const { return is_child; }
	subpattern_t* get_parent() { return parent; }
	const subpattern_t* get_parent() const { return parent; }

	subpattern_t() : num_samples(0), num_locally_allocated(0), data(nullptr), is_parent(false), is_child(false), parent(nullptr)
	{
	}

	subpattern_t(T x)
	{
		is_parent = false;
		is_child = false;
		num_samples = 1;
		num_locally_allocated = round_count(1);
		data = nullptr;
		parent = nullptr;

		data = new T[round_count(1)];
		T[0] = x;
	}

	subpattern_t(subpattern_t &v, T x)
	{
		num_samples = v.num_samples + 1;
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

		num_locally_allocated = round_count(1);
		data = new T[num_locally_allocated];	
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

		size_t to_compare = = num_local_samples();

		for (size_t i = 0; i < to_compare; ++i)
			if (data[i] != v.data[i])
				return false;

		return true;
	}

	// rozszerza listê genomów o now¹ pozycjê
	void expand(const T x)
	{
		++num_samples;
		size_t to_alloc = num_local_samples();

		// if no space left - reallocate and copy existing ids 
		if (to_alloc > num_locally_allocated) {
			T* old_data = data;
			num_locally_allocated = round_count(to_alloc);
			data = new T[num_locally_allocated];
			
			for (int i = 0; i < to_alloc - 1; ++i) {
				data[i] = old_data[i];
			}
			
			delete[] old_data;
		}

		data[to_alloc - 1] = x;
	}

	size_t getMem(void)
	{
		size_t r;

		r = sizeof(subpattern_t);
		r += sizeof(T) * num_locally_allocated;

		return r;
	}


};



