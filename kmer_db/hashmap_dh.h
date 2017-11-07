#pragma once

// Szablon klasy tablicy haszuj¹cej wykorzystujacy podwójne mieszanie

#include <mmintrin.h>
#include <cstdint>
#include <xmmintrin.h>
#include <iostream> 
#include <cstddef>

template<typename T>
inline size_t my_hasher_dh1(T x)
{
	return 0;			// !!! Fake impl.
}


template<>
inline size_t my_hasher_dh1<uint64_t>(uint64_t x)
{
	return x * 0xc70f6907ull;
}


template<>
inline size_t my_hasher_dh1<uint32_t>(uint32_t x)
{
	return x * 0xc70f6907ul;
}

template<typename T>
inline size_t my_hasher_dh2(T x)
{
	return 1;			// !!! Fake impl.
}


template<>
inline size_t my_hasher_dh2<uint64_t>(uint64_t x)
{
	return x * 0xc7096f07ull;
}


template<>
inline size_t my_hasher_dh2<uint32_t>(uint32_t x)
{
	return x * 0xc7096f07ull;
}


template <typename Key, typename Value>
class hash_map_dh {
public:
	typedef struct {
		Key key;
		Value val;
	} item_t;

private:
	Key empty_key;
	double max_fill_factor;

	size_t size;
	size_t filled;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;
	size_t allocated_mask2;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	void restruct(void)
	{
		item_t *old_data = data;
		size_t old_allocated = allocated;

		if (filled > old_allocated * max_fill_factor)
			allocated *= 2;

		allocated_mask = allocated - 1ull;
		allocated_mask2 = allocated_mask >> 1;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		clear();

		ht_memory += allocated * sizeof(item_t);

		//cout << "\n--- Realloc to: " << allocated << endl;

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].key != empty_key)
				insert(old_data[i].key, old_data[i].val);

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}



public:
	// fixme: iterator-like functionality - change to iterator
	item_t* begin() { return data; }
	const item_t* cbegin() const { return data; }
	item_t* end() { return data + allocated; }
	const item_t* cend() const { return data + allocated; }

	bool is_free(const item_t& item) const {
		return item.key == empty_key;
	}


	hash_map_dh(Key k_empty)
	{
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = 16;
		allocated_mask = allocated - 1;
		allocated_mask2 = allocated_mask >> 1;

		size = 0;
		filled = 0;
		data = new item_t[allocated];
		max_fill_factor = 0.8;
		
		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);

		set_special_keys(k_empty);
	}

	~hash_map_dh()
	{
		if (data)
			delete[] data;
	}

	size_t get_bytes() const {
		return ht_memory;
	}

	void set_special_keys(Key k_empty)
	{
		empty_key = k_empty;
		
		clear();
	}

	void clear(void)
	{
		size = 0;
		filled = 0;
		for (size_t i = 0; i < allocated; ++i)
			data[i].key = empty_key;
	}

	// Mozna to przyspieszyc tak, zebyinsert wykorzystywal wiedze o tym gdzie skonczyl szukac find
	Value* insert(Key k, Value v)
	{
		if (size >= size_when_restruct)
			restruct();

		size_t h = my_hasher_dh1<Key>(k) & allocated_mask;

		if (data[h].key != empty_key)
		{
			size_t h_step = 2ull * (my_hasher_dh2<Key>(k) & allocated_mask2) + 1ull;

			do
			{
				h += h_step;
				h = h & allocated_mask;
			} while (data[h].key != empty_key);
		}

		if (data[h].key == empty_key)
			++size;

		++filled;

		data[h].key = k;
		data[h].val = v;

		return &(data[h].val);
	}

	const Value* cfind(Key k) const {
		return find(k);
	}

	Value* find(Key k) const
	{
		size_t h = my_hasher_dh1<Key>(k) & allocated_mask;
		if (data[h].key == k)
			return &(data[h].val);

		if (data[h].key == empty_key)
			return nullptr;

		size_t h_step = 2ull * (my_hasher_dh2<Key>(k) & allocated_mask2) + 1ull;
		h += h_step;
		h = h & allocated_mask;

		//		++ht_total;
		while (data[h].key != empty_key)
		{
			if (data[h].key == k)
			{
				//				++ht_match;
				return &(data[h].val);
			}
			else
			{
				h += h_step;
				h = h & allocated_mask;
			}
		}

		return nullptr;
	}

	//	__declspec(noinline) 
	void prefetch(Key k)
	{
		size_t h = my_hasher_dh1<Key>(k) & allocated_mask;
		size_t h_step = 2ull * (my_hasher_dh2<Key>(k) & allocated_mask2) + 1ull;

#ifdef WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
		h = (h + h_step) & allocated_mask;
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
		h = (h + h_step) & allocated_mask;
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const
	{
		return filled;
	}

	void reserve_for_additional(size_t n_elems)
	{
		if (filled + n_elems <= allocated * max_fill_factor)
			return;


		item_t *old_data = data;
		size_t old_allocated = allocated;

		while (filled + n_elems > allocated * max_fill_factor)
			allocated *= 2;

		allocated_mask = allocated - 1ull;
		allocated_mask2 = allocated_mask >> 1;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		std::cout << "\n--- Realloc to: " << allocated << "...";

		data = new item_t[allocated];
		clear();

		ht_memory += allocated * sizeof(item_t);

		std::cout << "done!" << std::endl;

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].key != empty_key)
				insert(old_data[i].key, old_data[i].val);

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

};