#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

// *****************************************************************************************
//
template<typename T, unsigned ALIGNMENT = 4096> class aligned_vector
{
	T *data_aligned;
	char *data_raw;
	size_t data_size;
	size_t data_allocated;

	// *****************************************************************************************
	//
	void allocate(void)
	{
		if (data_raw)
			delete[] data_raw;

		data_allocated = data_size;
		if (data_allocated < 16)
			data_allocated = 16;

		size_t bytes_to_allocate = data_allocated * sizeof(T) + ALIGNMENT;
		
		data_raw = new char[bytes_to_allocate];

		size_t address = reinterpret_cast<size_t>(data_raw);
		size_t offset = ALIGNMENT - address % ALIGNMENT;

		data_aligned = reinterpret_cast<T*>(data_raw + offset);
	}

	// *****************************************************************************************
	//
	void free()
	{
		if (data_raw)
			delete[] data_raw;
		data_raw = nullptr;
	}

public:
	typedef T value_type;

	size_t get_bytes() const {
		return data_allocated * sizeof(T) + ALIGNMENT;
	}

	// *****************************************************************************************
	//
	aligned_vector(size_t _data_size = 0) : data_size(_data_size), data_aligned(nullptr), data_raw(nullptr)
	{
		allocate();
	}

	// *****************************************************************************************
	//
	~aligned_vector()
	{
		free();
	}

	// *****************************************************************************************
	//
	void swap(aligned_vector<T, ALIGNMENT> &x)
	{
		::swap(data_aligned, x.data_aligned);
		::swap(data_raw, x.data_raw);
		::swap(data_size, x.data_size);
		::swap(data_allocated, x.data_allocated);
	}

	// *****************************************************************************************
	//
	size_t size()
	{
		return data_size;
	}

	// *****************************************************************************************
	//
	void resize(size_t new_size)
	{
		if (new_size > data_allocated)
		{
			data_size = new_size;
			allocate();
		}
		else
			data_size = new_size;
	}

	// *****************************************************************************************
	//
	T* begin()		// pseudo iterator
	{
		return data_aligned;
	}

	// *****************************************************************************************
	//
	T* end()		// pseudo iterator
	{
		return data_aligned + data_size;
	}

	// *****************************************************************************************
	//
	T* data()
	{
		return data_aligned;
	}

	// *****************************************************************************************
	//
	T& operator[](size_t pos)
	{
		return data_aligned[pos];
	}

	// *****************************************************************************************
	//
	const T& operator[](size_t pos) const
	{
		return data_aligned[pos];
	}
};