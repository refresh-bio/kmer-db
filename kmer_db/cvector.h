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
