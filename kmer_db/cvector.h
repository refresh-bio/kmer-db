#pragma once
#include <algorithm>

using namespace std;

size_t cv_memory = 0;

size_t mem_pattern_desc = 0;		// iloœæ pamiêci zajmowana przez wszystkie wzorce
size_t mem_os_pattern_desc = 0;		// j.w., ale uwzglêdnia to, ¿e OS bêdzie alokowa³ pamiêæ dynamiczn¹ z wyr. do 16B

template<typename T>
class const_vector {
	size_t size;
	T* data;

public:
	const_vector()
	{
		size = 0;
		data = nullptr;
	}

	const_vector(T x)
	{
		size = 1;
		data = new T[1];
		data[0] = x;

		cv_memory += sizeof(T);

//		++cv_objs;
	}

	const_vector(const const_vector &v, T x)
	{
		size = v.size + 1;
		data = new T[size];
//		++cv_objs;
		cv_memory += sizeof(T) * size;

		if(size > 1)
			copy_n(v.data, size - 1, data);
		data[size - 1] = x;
	}

	const_vector(const const_vector &v)
	{
		size = v.size;
		if (size)
		{
			data = new T[size];
			copy_n(v.data, size, data);

//			++cv_objs;
			cv_memory += sizeof(T) * size;
		}
		else
			data = nullptr;
	}

	const_vector(const_vector &&v)
	{
		size = v.size;
		data = v.data;

		v.size = 0;
		v.data = nullptr;
	}

	~const_vector()
	{
		if (data)
		{
			delete[] data;
//			--cv_objs;
			cv_memory -= sizeof(T) * size;
		}
	}

	const_vector<T>& operator=(const const_vector<T> &v)
	{
		if (this == &v)
			return *this;

		if (data)
		{
			delete[] data;
//			--cv_objs;
			cv_memory -= sizeof(T) * size;
		}

		size = v.size;
		data = new T[size];
		copy_n(v.data, size, data);

//		++cv_objs;
		cv_memory += sizeof(T) * size;

		return *this;
	}

	const_vector<T>& operator=(const_vector<T> &&v)
	{
		if (this == &v)
			return *this;

		if (data)
		{
			delete[] data;
//			--cv_objs;
			cv_memory -= sizeof(T) * size;
		}

		size = v.size;
		data = v.data;

		v.size = 0;
		v.data = nullptr;

		return *this;
	}

	bool operator==(const const_vector<T> &v)
	{
		if (v.size != size)
			return false;

		for (size_t i = 0; i < size; ++i)
			if (data[i] != v.data[i])
				return false;

		return true;
	}

	size_t get_size(void)
	{
		return size;
	}
};


template<typename T>
class fake_const_vector {
	size_t size;
	T* data;

public:
	fake_const_vector()
	{
		size = 0;
		data = nullptr;
	}

	fake_const_vector(T x)
	{
		size = 1;
		data = nullptr;

		cv_memory += sizeof(fake_const_vector);

		//		++cv_objs;
	}

	fake_const_vector(const fake_const_vector &v, T x)
	{
		size = v.size + 1;
		data = v.data;
		//		++cv_objs;
		cv_memory += sizeof(fake_const_vector);
	}

	fake_const_vector(const fake_const_vector &v)
	{
		size = v.size;
		if (size)
		{
			data = v.data;
			cv_memory += sizeof(fake_const_vector);
		}
		else
			data = nullptr;
	}

	fake_const_vector(fake_const_vector &&v)
	{
		size = v.size;
		data = v.data;

		v.size = 0;
		v.data = nullptr;
	}

	~fake_const_vector()
	{
		cv_memory -= sizeof(fake_const_vector);
	}

	fake_const_vector<T>& operator=(const fake_const_vector<T> &v)
	{
		if (this == &v)
			return *this;

		cv_memory += sizeof(fake_const_vector);

		return *this;
	}

	fake_const_vector<T>& operator=(fake_const_vector<T> &&v)
	{
		if (this == &v)
			return *this;

		size = v.size;
		data = v.data;

		v.size = 0;
		v.data = nullptr;

		return *this;
	}

	bool operator==(const fake_const_vector<T> &v)
	{
		if (v.size != size)
			return false;

/*		for (size_t i = 0; i < size; ++i)
			if (data[i] != v.data[i])
				return false;*/

		return true;
	}

	size_t get_size(void)
	{
		return size;
	}
};


template<typename T>
class pattern_desc {
	bool is_parrent;			// informacja czy wêze³ jest rodzicem innego, tzn. czy jakiœ wêze³ jest zapisany jako rozszerzenie bie¿¹cego [true]
	bool is_child;				// informacja czy wêze³ jest zapisany jako rozszerzenie innego [true] czy te¿ wprost (lista id genomów) [false]
	size_t no_ids;				// liczba id genomów w wêŸle	
	T* data;					// tablica id genomów
	pattern_desc* parrent;		// wskaŸnik do wêz³a rodzica

	size_t round_up_16(size_t x)
	{
		if (x % 16)
			x = ((x / 16) + 1) * 16;

		return x;
	}

public:
	pattern_desc() : no_ids(0), data(nullptr), is_parrent(false), is_child(false), parrent(nullptr)
	{
		int aa = 1;

	}

	pattern_desc(T x)
	{
		is_parrent = false;
		is_child = false;
		no_ids = 1;
		data = nullptr;
		parrent = nullptr;

		mem_pattern_desc += sizeof(T);
		mem_os_pattern_desc += round_up_16(sizeof(T));
	}

	pattern_desc(pattern_desc &v, T x)
	{
		no_ids = v.no_ids + 1;
		data = nullptr;

		if (v.no_ids == 0)
		{
			is_parrent = false;
			is_child = false;
			v.is_parrent = false;
		}
		else
		{
			is_parrent = false;
			is_child = true;
			v.is_parrent = true;
			parrent = &v;
		}

		size_t to_alloc = 1;

		mem_pattern_desc += sizeof(T);
		mem_os_pattern_desc += round_up_16(sizeof(T));
	}

	pattern_desc(const pattern_desc &v)
	{
		is_parrent = v.is_parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		data = nullptr;
		parrent = v.parrent;

		if (no_ids)
		{
			size_t to_alloc = no_ids;
			
			if(parrent)
				to_alloc -= parrent->no_ids;

			mem_pattern_desc += sizeof(T) * to_alloc;
			mem_os_pattern_desc += round_up_16(sizeof(T) * to_alloc);
		}
		else
			data = nullptr;
	}

	pattern_desc(pattern_desc &&v)
	{
		is_parrent = v.is_parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		data = v.data;
		parrent = v.parrent;

		v.data = nullptr;
		v.parrent = nullptr;
		v.no_ids = 0;
	}

	~pattern_desc()
	{
		if (no_ids)
		{
//			delete[] data;

			size_t to_free = no_ids;
			if (parrent)
				to_free -= parrent->no_ids;

			mem_pattern_desc -= sizeof(T) * to_free;
			mem_os_pattern_desc -= round_up_16(sizeof(T) * to_free);
		}
	}

	pattern_desc<T>& operator=(const pattern_desc<T> &v)
	{
		if (this == &v)
			return *this;

		if (no_ids)
		{
//			delete[] data;

			size_t to_free = no_ids;
			if (parrent)
				to_free -= parrent->no_ids;

			mem_pattern_desc -= sizeof(T) * to_free;
			mem_os_pattern_desc -= round_up_16(sizeof(T) * to_free);
		}

		is_parrent = v.parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		parrent = v.parrent;
		data = nullptr;

		size_t to_alloc = no_ids;

		if (parrent)
			to_alloc -= parrent->no_ids;

		mem_pattern_desc += sizeof(T) * to_alloc;
		mem_os_pattern_desc += round_up_16(sizeof(T) * to_alloc);

		return *this;
	}

	pattern_desc<T>& operator=(pattern_desc<T> &&v)
	{
		if (this == &v)
			return *this;

		if (no_ids)
		{
			//			delete[] data;

			size_t to_free = no_ids;
			if (parrent)
				to_free -= parrent->no_ids;

			mem_pattern_desc -= sizeof(T) * to_free;
			mem_os_pattern_desc -= round_up_16(sizeof(T) * to_free);
		}

		is_parrent = v.parrent;
		is_child = v.is_child;
		no_ids = v.no_ids;
		parrent = v.parrent;
		data = nullptr;

		size_t to_alloc = no_ids;

		if (parrent)
			to_alloc -= parrent->no_ids;

		mem_pattern_desc += sizeof(T) * to_alloc;
		mem_os_pattern_desc += round_up_16(sizeof(T) * to_alloc);

		v.data = nullptr;
		v.parrent = nullptr;
		v.no_ids = 0;

		return *this;
	}

	bool operator==(const pattern_desc<T> &v)
	{
		if (v.no_ids != no_ids)
			return false;
		if (v.parrent != parrent)
			return false;

		size_t to_compare = no_ids;
		if (parrent)
			to_compare -= parrent->no_ids;

		for (size_t i = 0; i < to_compare; ++i)
			if (data[i] != v.data[i])
				return false;

		return true;
	}

	// rozszerza listê genomów o now¹ pozycjê
	void expand(const T x)
	{
		++no_ids;
		size_t to_alloc = no_ids;
		if (parrent)
			to_alloc -= parrent->no_ids;
		
		T* old_data = data;
		data = nullptr;
//		data = new T[to_alloc];
		// kopiowanie

		mem_pattern_desc += sizeof(T);
		mem_os_pattern_desc += round_up_16(sizeof(T) * no_ids);
		mem_os_pattern_desc -= round_up_16(sizeof(T) * (no_ids-1));
	}

	size_t size(void)
	{
		size_t r;

		r = sizeof(pattern_desc);
		r += sizeof(T) * no_ids;
		if (parrent)
			r -= sizeof(T) * parrent->no_ids;

		return r;
	}

	bool get_is_parrent()
	{
		return is_parrent;
	}

	bool get_is_child()
	{
		return is_child;
	}
};


