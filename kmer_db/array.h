#pragma once
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

template <class T>
class Array
{

public:
	size_t size() const { return width; }
	const std::vector<T>& getData() const { return v; }
	std::vector<T>& getData() { return v; }

	Array(int size) : width(size), height(size), v(size * size) {}
	Array() : width(0), height(0) {}
	Array(const Array<T>& ref) : width(ref.width), height(ref.height), v(ref.v) {}
	Array(Array&& rhs) : width(rhs.width), height(rhs.height), v(std::move(rhs.v)) {}

	Array<T>& operator=(const Array<T>& ref) = default;
	Array<T> operator-(const Array<T>& b) {
		Array<T> out;
		if (this->size() == b.size()) {
			out.resize(this->size(), this->size(), 0);

			for (int i = 0; i < this->size() * this->size(); ++i) {
				out.getData()[i] = this->getData()[i] - b.getData()[i];
			}
		}
		return out;
	}

	void resize(int width, int height) {
		this->width = width;
		this->height = height;
		v.resize(width * height);
	}

	void resize(int width, int height, const T& value) {
		this->width = width;
		this->height = height;
		v.resize(width * height, value);
	}

	void clear() {
		v.clear();
	}

	T* operator[](const int row) { return v.data() + row * width; }
	const T* operator[](const int row) const { return v.data() + row * width; }

	void save(std::ofstream & file) {
		for (int i = 0; i < size(); ++i) {
			for (int j = 0; j < size(); ++j) {
				file  << (*this)[i][j] << ',';
			}
			file << std::endl;
		}
	}


protected:
	int width;
	int height;
	std::vector<T> v;
};

template <class T>
class LowerTriangularMatrix {
public:
	size_t getSize() const { return size; }
	
	LowerTriangularMatrix() : size(0) {}
	LowerTriangularMatrix(size_t size) : size(size), data(size * (size - 1) / 2) {}

	T* operator[](size_t i) { return data.data() + i * (i - 1) / 2; }
	const T* operator[](size_t i) const { return data.data() + i * (i - 1) / 2; }

	T* at(size_t i, size_t j) { return data[i * (i - 1) / 2  + j]; }
	const T* at(size_t i, size_t j) const { return data[i * (i - 1) / 2 + j]; }

	void resize(size_t size) {
		this->size = size;
		data.resize(size * (size - 1) / 2);
	}

	void clear() {
		data.clear();
	}

	void save(std::ofstream & file) {
		T * ptr = data.data();
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < i; ++j) {
				file << *ptr++ << ',';
			}

			for (int j = i; j < size; ++j) {
				file << "0,";
			}

			file << std::endl;
		}
	}

protected:
	size_t size;
	std::vector<T> data;

};


class ArrayBuffer {
	int32_t *buf;
	int32_t *cur_buf;
	int32_t *end_buf;
	uint32_t *vec;
	size_t max_buf_size;
	size_t buf_size;

//	__declspec(noinline)
	void update()
	{
		int32_t cnt = 0;

/*		for (size_t i = 0; i < buf_size; ++i)
			if (buf[i] < 0)
				cnt = -buf[i];
			else
				vec[buf[i]] += cnt;

		buf[0] = -cnt;
		buf_size = 1;*/

		for (int32_t *p = buf; p != cur_buf; ++p)
			if (*p < 0)
				cnt = -*p;
			else
				vec[*p] += cnt;

		cur_buf = buf;
		*cur_buf++ = -cnt;
	}

public:
	ArrayBuffer() : vec(nullptr), buf(nullptr)
	{};

	~ArrayBuffer()
	{
		if(buf)
			delete[] buf;
	}

	void Assign(uint32_t *_vec, size_t _max_buf_size)
	{
		vec = _vec;
		max_buf_size = _max_buf_size;
		buf_size = 0;
		buf = new int32_t[max_buf_size];
		cur_buf = buf;
		end_buf = buf + max_buf_size;
	}

	//__declspec(noinline)
	void Push(uint32_t id)
	{
//		if (buf_size >= max_buf_size)
		if(cur_buf == end_buf)
			update();

//		buf[buf_size++] = (int32_t)id;
//		*cur_buf++ = (int32_t)id;
		*cur_buf++ = id;
	}

//	__declspec(noinline)
	void SetCounter(uint32_t cnt)
	{
		if (cur_buf == end_buf)
//		if (buf_size >= max_buf_size)
			update();
		
//		buf[buf_size++] = -(int32_t) cnt;
		*cur_buf++ = -(int32_t)cnt;
	}

	void Prefetch()
	{
#ifdef WIN32
//		_mm_prefetch((const char*)(buf + buf_size), _MM_HINT_T0);
		_mm_prefetch((const char*) cur_buf, _MM_HINT_T0);
#else
//		__builtin_prefetch(buf + buf_size);
		__builtin_prefetch(cur_buf);
#endif
	}

	void Finish()
	{
		update();
	}
};
