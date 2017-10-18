#pragma once

// Szablon klasy tablicy haszuj¹cej wykorzystujacy adresowanie liniowe

#include <mmintrin.h>

size_t ht_memory = 0;

size_t ht_total = 0;
size_t ht_match = 0;

template<typename T>
size_t my_hasher(T x)
{
	return 0;			// !!! Fake impl.
}

template<>
size_t my_hasher<uint64_t>(uint64_t x)
{
	return x * 0xc70f6907ull;
}

template<>
size_t my_hasher<uint32_t>(uint32_t x)
{
	return x * 0xc70f6907ul;
}

template <typename Key, typename Value>
class hash_map {
	typedef struct {
		Key key;
		Value val;
	} item_t;
	
	Key empty_key;
	Key erased_key;
	double max_fill_factor;
	double max_rest_factor;

	size_t size;
	size_t filled;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	void restruct(void)
	{
		item_t *old_data = data;
		size_t old_allocated = allocated;

		if(filled > old_allocated * max_rest_factor)
			allocated *= 2;

		allocated_mask = allocated - 1;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		clear();

		ht_memory += allocated * sizeof(item_t);

		cout << "\n--- Realloc to: " << allocated << endl;

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].key != empty_key && old_data[i].key != erased_key)
				insert(old_data[i].key, old_data[i].val);

		delete[] old_data;
		ht_memory -= old_allocated  * sizeof(item_t);
	}

public:
	hash_map()
	{
		allocated = 16;
		allocated_mask = allocated - 1;

		size = 0;
		filled = 0;
		data = new item_t[allocated];
		max_fill_factor = 0.45;
		max_rest_factor = 0.35;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);
	}

	~hash_map()
	{
		if (data)
			delete[] data;
	}

	void set_special_keys(Key k_empty, Key k_erased)
	{
		empty_key = k_empty;
		erased_key = k_erased;

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

		size_t h = my_hasher<Key>(k) & allocated_mask;

		while (data[h].key != empty_key && data[h].key != erased_key)
			h = (h + 1) & allocated_mask;

		if (data[h].key == empty_key)
			++size;

		++filled;

		data[h].key = k;
		data[h].val = v;

		return &(data[h].val);
	}

	Value* find(Key k)
	{
		size_t h = my_hasher<Key>(k) & allocated_mask;

//		++ht_total;
		while (data[h].key != empty_key)
		{
			if (data[h].key == k)
			{
//				++ht_match;
				return &(data[h].val);
			}
			else
				if (h == allocated_mask)
					h = 0;
				else
					++h;
			//				h = (h + 1) & allocated_mask;

//			++ht_total;
		}

		return nullptr;
	}

//	__declspec(noinline) 
	void prefetch(Key k)
	{
		size_t h = my_hasher<Key>(k) & allocated_mask;

#ifdef WIN32
		_mm_prefetch((const char*) (data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}
		
	void erase(Key k)
	{
		size_t h = my_hasher<Key>(k) & allocated_mask;

		while (data[h].key != empty_key)
			if (data[h].key == k)
			{
				data[h].key = erased_key;
				data[h].val = Value();
				--filled;
				return;
			}
			else
				h = (h + 1) & allocated_mask;
	}

	size_t get_size(void)
	{
		return filled;
	}

	void reserve_for_additional(size_t n_elems)
	{
		if (filled + n_elems <= allocated * max_rest_factor)
			return;


		item_t *old_data = data;
		size_t old_allocated = allocated;

		while (filled + n_elems > allocated * max_rest_factor)
			allocated *= 2;

		allocated_mask = allocated - 1;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		cout << "\n--- Realloc to: " << allocated << "...";

		data = new item_t[allocated];
		clear();

		ht_memory += allocated * sizeof(item_t);

		cout << "done!" << endl;

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].key != empty_key && old_data[i].key != erased_key)
				insert(old_data[i].key, old_data[i].val);

		delete[] old_data;
		ht_memory -= old_allocated  * sizeof(item_t);
	}

};