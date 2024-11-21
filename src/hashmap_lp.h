#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#ifdef ARCH_X64
#include <mmintrin.h>
#include <xmmintrin.h>
#else
#include <arm_neon.h>
#endif

#include <cstdint>
#include <iostream> 
#include <cstddef>
#include <thread>
#include <fstream>

#include "log.h"

// *****************************************************************************************
//
template<typename T>
inline size_t my_hasher_lp(T x)
{
	return 0;			// !!! Fake impl.
}

// *****************************************************************************************
//
template<>
inline size_t my_hasher_lp<uint64_t>(uint64_t h)
{
//	return x * 0xc70f6907ull;

	h ^= h >> 33;
	h *= 0xff51afd7ed558ccdL;
	h ^= h >> 33;
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= h >> 33;

	return h;

}

// *****************************************************************************************
//
template<>
inline size_t my_hasher_lp<uint32_t>(uint32_t h)
{
//	return h * 0xc70f6907ul;

	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return h;
}

// *****************************************************************************************
//
template <typename Key, typename Value>
class hash_map_lp {
public:
	typedef struct {
		Key key;
		Value val;
	} item_t;

	static const size_t INITIAL_SIZE = 16;
//	static const Key empty_key = static_cast<Key>(-1);
	static const Value empty_value = std::numeric_limits<Value>::max();	

private:
	double max_fill_factor;

	size_t filled;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	
	// *****************************************************************************************
	//
	void restruct(void)
	{
		item_t *old_data = data;
		size_t old_allocated = allocated;

		if (filled >= size_when_restruct)
			allocated *= 2;

		allocated_mask = allocated - 1ull;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		clear();

		ht_memory += allocated * sizeof(item_t);
		filled = 0;

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].val != empty_value)
				fast_insert(old_data[i].key, old_data[i].val);

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

public:
	// *****************************************************************************************
	//
	// fixme: iterator-like functionality - change to iterator
	item_t* begin() { return data; }
	const item_t* cbegin() const { return data; }
	item_t* end() { return data + allocated; }
	const item_t* cend() const { return data + allocated; }

	
	size_t get_size(void) const { return filled; }
	size_t get_capacity(void) const { return allocated; }
	bool is_free(const item_t& item) const { return item.val == empty_value; }
	


	// *****************************************************************************************
	//
	hash_map_lp()
	{
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = INITIAL_SIZE;
		allocated_mask = allocated - 1;
	
		filled = 0;
		data = new item_t[allocated];
		max_fill_factor = 0.8;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);

		clear();
	}

	// *****************************************************************************************
	//
	~hash_map_lp()
	{
		if (data)
			delete[] data;
	}

	// *****************************************************************************************
	//
	size_t get_bytes() const {
		return ht_memory;
	}


	// *****************************************************************************************
	//
	void clear(void)
	{
		filled = 0;
		for (size_t i = 0; i < allocated; ++i)
		{
			data[i].key = Key{};
			data[i].val = empty_value;
		}
	}

	// *****************************************************************************************
	//
	void release()
	{
		if (data)
			delete[] data;

		data = nullptr;
	}

	// *****************************************************************************************
	//
	void parallel_clear(void)
	{
		filled = 0;

		int n_threads = std::max((int)std::thread::hardware_concurrency(), 1);
		std::vector<std::thread> threads(n_threads);

		//LOG_DEBUG << "Clearing hashtable (parallel)..." ;

		for (int tid = 0; tid < n_threads; ++tid) {
			threads[tid] = std::thread([tid, n_threads, this] {
				size_t block = allocated / n_threads;
				size_t lo = tid * block;
				size_t hi = (tid == n_threads - 1) ? allocated : lo + block;

				item_t * p = data + lo;
				item_t * end = data + hi;

				for (; p < end; ++p) {
					p->key = Key{};
					p->val = empty_value;
				}

			});
		}

		for (auto& t : threads) {
			t.join();
		}

	}

	// *****************************************************************************************
	//
	uint32_t hash_val32(const Key k)
	{
		return (uint32_t) (my_hasher_lp<Key>(k) & allocated_mask);
	}

	// *****************************************************************************************
	//
	uint64_t hash_val64(const Key k)
	{
		return (uint64_t) (my_hasher_lp<Key>(k) & allocated_mask);
	}

	// *****************************************************************************************
	//
	Value* insert(Key k, Value v)
	{
		if (filled >= size_when_restruct) {
			throw std::runtime_error("Assertion error : hashmap_lp::restruct() should never be invoked");
			//restruct();
		}

		size_t h = my_hasher_lp<Key>(k) & allocated_mask;

		if (data[h].val != empty_value)
		{
			do
			{
				h = (h + 1) & allocated_mask;
			} while (data[h].val != empty_value);
		}

		++filled;

		data[h].key = k;
		data[h].val = v;

		return &(data[h].val);
	}

	// *****************************************************************************************
	//
	void fast_insert(Key k, Value v)
	{
		size_t h = my_hasher_lp<Key>(k) & allocated_mask;

		if (data[h].val != empty_value)
		{
			do
			{
				h = (h + 1) & allocated_mask;
			} while (data[h].val != empty_value);
		}

		++filled;

		data[h].key = k;
		data[h].val = v;
	}

	// *****************************************************************************************
	//
	void insert(Key k, Value v, item_t* place) 
	{
		place->key = k;
		place->val = v;
		++filled;
	}

	// *****************************************************************************************
	//
	const Value* cfind(Key k) const {
		return find(k);
	}

	// *****************************************************************************************
	//
	Value* find(Key k) const
	{
		size_t h = my_hasher_lp<Key>(k) & allocated_mask;

		if (data[h].val == empty_value)
			return nullptr;

		if (data[h].key == k)
			return &(data[h].val);

		h = (h + 1) & allocated_mask;

		while (data[h].val != empty_value)
		{
			if (data[h].key == k)
			{
				return &(data[h].val);
			}
			else
			{
				h = (h + 1) & allocated_mask;
			}
		}

		return nullptr;
	}


	// *****************************************************************************************
	//
	item_t* find_item(Key k) const
	{
		size_t h = my_hasher_lp<Key>(k) & allocated_mask;
		
//		while (data[h].key != k && data[h].val != empty_value) {
		while (data[h].val != empty_value && data[h].key != k) {
			h = (h + 1) & allocated_mask;
		} 
		return &(data[h]);

/*		if (data[h].val == empty_value || data[h].key == k)
			return &(data[h]);

		for (++h; h < allocated; ++h)
			if(data[h].val == empty_value || data[h].key == k)
				return &(data[h]);

		for (h = 0; h < allocated; ++h)
			if (data[h].val == empty_value || data[h].key == k)
				return &(data[h]);

		return &(data[h]);*/
	}

	// *****************************************************************************************
	//
	void prefetch(Key k) const
	{
		size_t h = my_hasher_lp<Key>(k) & allocated_mask;

#ifdef WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	// *****************************************************************************************
	//
	Value& operator[](const Key key)
	{
		if (filled >= size_when_restruct) {
			restruct();
		}

		size_t h = my_hasher_lp<Key>(key) & allocated_mask;

		if (data[h].val == empty_value)
		{
			data[h].key = key;
			data[h].val = 0;
			++filled;

			return data[h].val;
		}

		if (data[h].key == key)
			return data[h].val;

		h = (h + 1) & allocated_mask;

		while (data[h].val != empty_value)
		{
			if (data[h].key == key)
			{
				return data[h].val;
			}
			else
			{
				h = (h + 1) & allocated_mask;
			}
		}

		data[h].key = key;
		data[h].val = 0;
		++filled;

		return data[h].val;
	}

	// *****************************************************************************************
	//
	bool need_reserve_for_additional(size_t n_elems)
	{
		return filled + n_elems > allocated * max_fill_factor;
	}

	// *****************************************************************************************
	//
	bool reserve_for_additional(size_t n_elems)
	{
		if (filled + n_elems <= allocated * max_fill_factor)
			return false;

		item_t *old_data = data;
		size_t old_allocated = allocated;

	//	LOG_DEBUG << "reserve_for_additional - in\n" ;
		while (filled + n_elems > allocated * max_fill_factor)
			allocated *= 2;

	//	LOG_DEBUG << "reserve_for_additional - new_size: " << allocated << std::endl ;

		allocated_mask = allocated - 1ull;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

	//	LOG_DEBUG << "\n--- Realloc to: " << allocated << "..." << std::endl ;

		data = new item_t[allocated];
	//	LOG_DEBUG << "reserve_for_additional - after new: " << allocated << std::endl ;

		clear();
	//	LOG_DEBUG << "reserve_for_additional - after clear: " << allocated << std::endl ;

		ht_memory += allocated * sizeof(item_t);

		for (size_t i = 0; i < old_allocated; ++i)
		{
			if (old_data[i].val != empty_value)
				insert(old_data[i].key, old_data[i].val);
		}

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);

		return true;
	}

	// *****************************************************************************************
	//
	void serialize(std::ofstream& file) const {
		
		size_t buf_size = (1 << 25) / sizeof(item_t);
		item_t* buffer = new item_t[buf_size];
		
		serialize(file, buffer, buf_size);

		delete [] buffer;
	}


	// *****************************************************************************************
	//
	void serialize(std::ofstream& file, item_t* buffer, size_t buf_size) const {
		
		file.write(reinterpret_cast<const char*>(&max_fill_factor), sizeof(max_fill_factor));

		file.write(reinterpret_cast<const char*>(&filled), sizeof(filled));
		file.write(reinterpret_cast<const char*>(&allocated), sizeof(allocated));
		file.write(reinterpret_cast<const char*>(&size_when_restruct), sizeof(size_when_restruct));
		file.write(reinterpret_cast<const char*>(&allocated_mask), sizeof(allocated_mask));
	
		file.write(reinterpret_cast<const char*>(&ht_memory), sizeof(ht_memory));
		file.write(reinterpret_cast<const char*>(&ht_total), sizeof(ht_total));
		file.write(reinterpret_cast<const char*>(&ht_match), sizeof(ht_match));

		// write bit vectors of filled slots
		size_t bv_size = (allocated + 63) / 64;
		uint64_t* bv = new uint64_t[bv_size];	
		std::fill_n(bv, bv_size, 0);

		for (uint64_t i = 0; i < allocated; ++i) {
			if (data[i].val != empty_value) {
				uint64_t word = i >> 6;
				uint64_t offset = i & 63;

				bv[word] |= 1ull << offset;
			}
		}
		file.write(reinterpret_cast<const char*>(bv), sizeof(uint64_t) * bv_size);

		// write HT filled slots in blocks
		size_t out_id = 0;

		for (uint64_t i = 0; i < allocated; ++i) {
			if (data[i].val != empty_value) {
				buffer[out_id++] = data[i];
			}

			// if buffer full
			if (out_id == buf_size) {
				file.write(reinterpret_cast<const char*>(buffer), sizeof(item_t) * buf_size);
				out_id = 0;
			}	
		}

		// write rest of the buffer
		file.write(reinterpret_cast<const char*>(buffer), sizeof(item_t) * out_id);
		
		delete[] bv;
	}

	// *****************************************************************************************
	//
	bool deserialize(std::ifstream& file) {

		size_t buf_size = (1 << 25) / sizeof(item_t);
		item_t* buffer = new item_t[buf_size];

		bool ok = deserialize(file, buffer, buf_size);

		delete[] buffer;

		return ok;
	}

	// *****************************************************************************************
	//
	bool deserialize(std::ifstream& file, item_t* buffer, size_t buf_size, bool skipData) {

		file.read(reinterpret_cast<char*>(&max_fill_factor), sizeof(max_fill_factor));

		file.read(reinterpret_cast<char*>(&filled), sizeof(filled));
		file.read(reinterpret_cast<char*>(&allocated), sizeof(allocated));
		file.read(reinterpret_cast<char*>(&size_when_restruct), sizeof(size_when_restruct));
		file.read(reinterpret_cast<char*>(&allocated_mask), sizeof(allocated_mask));
		
		file.read(reinterpret_cast<char*>(&ht_memory), sizeof(ht_memory));
		file.read(reinterpret_cast<char*>(&ht_total), sizeof(ht_total));
		file.read(reinterpret_cast<char*>(&ht_match), sizeof(ht_match));

		// load bit vectors of filled slots
		size_t bv_size = (allocated + 63) / 64;
		
		if (skipData) {
			// skip bit vectors and raw data
			file.seekg(sizeof(uint64_t) * bv_size, ios_base::cur);
			file.seekg(sizeof(uint64_t) * filled, ios_base::cur);
		}
		else {
			uint64_t* bv = new uint64_t[bv_size];

			file.read(reinterpret_cast<char*>(bv), sizeof(uint64_t) * bv_size);

			// free already allocated memory
			delete[] data;
			data = new item_t[allocated];

			// load in portions
			size_t buf_id = 0;
			size_t remaining = filled;

			for (uint64_t i = 0; i < allocated; ++i) {
				// load data 
				if (buf_id == 0 || buf_id == buf_size) {
					file.read(reinterpret_cast<char*>(buffer), sizeof(uint64_t) * std::min(buf_size, remaining));
					buf_id = 0;
					remaining -= std::min(buf_size, remaining);
				}

				uint64_t word = i >> 6;
				uint64_t offset = i & 63;

				// if slot is filled
				if (bv[word] & (1ull << offset)) {
					data[i] = buffer[buf_id++];
				}
				else {
					data[i].key = Key{};
					data[i].val = empty_value;
				}
			}

			delete[] bv;
		}

		return file.good();
	}

	// *****************************************************************************************
	//
	bool deserialize_into_vector(std::ifstream& file, item_t* buffer, size_t buf_size, vector<pair<Key, Value>> &vec, bool skipData) {

		file.read(reinterpret_cast<char*>(&max_fill_factor), sizeof(max_fill_factor));

		file.read(reinterpret_cast<char*>(&filled), sizeof(filled));
		file.read(reinterpret_cast<char*>(&allocated), sizeof(allocated));
		file.read(reinterpret_cast<char*>(&size_when_restruct), sizeof(size_when_restruct));
		file.read(reinterpret_cast<char*>(&allocated_mask), sizeof(allocated_mask));

		file.read(reinterpret_cast<char*>(&ht_memory), sizeof(ht_memory));
		file.read(reinterpret_cast<char*>(&ht_total), sizeof(ht_total));
		file.read(reinterpret_cast<char*>(&ht_match), sizeof(ht_match));

		// load bit vectors of filled slots
		size_t bv_size = (allocated + 63) / 64;

		if (skipData) {
			// skip bit vectors and raw data
			file.seekg(sizeof(uint64_t) * bv_size, ios_base::cur);
			file.seekg(sizeof(uint64_t) * filled, ios_base::cur);
		}
		else {
			file.seekg(sizeof(uint64_t) * bv_size, ios_base::cur);

			vec.resize(filled);

			// load in portions
			size_t remaining = filled;

			char* vec_ptr = reinterpret_cast<char*>(vec.data());

			while (remaining)
			{
				file.read(vec_ptr, sizeof(uint64_t) * std::min(buf_size, remaining));
				vec_ptr += sizeof(uint64_t) * std::min(buf_size, remaining);
				remaining -= std::min(buf_size, remaining);
			}
		}

		return file.good();
	}
};