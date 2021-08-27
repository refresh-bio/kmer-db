#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "conversion.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

using namespace std;

// *****************************************************************************************
//
template <class T>
class Array
{

public:
	size_t size() const { return width; }
	const std::vector<T>& getData() const { return v; }
	std::vector<T>& getData() { return v; }

	// *****************************************************************************************
	//
	Array(int size) : width(size), height(size), v(size * size) {}

	// *****************************************************************************************
	//
	Array() : width(0), height(0) {}

	// *****************************************************************************************
	//
	Array(const Array<T>& ref) : width(ref.width), height(ref.height), v(ref.v) {}

	// *****************************************************************************************
	//
	Array(Array&& rhs) : width(rhs.width), height(rhs.height), v(std::move(rhs.v)) {}

	// *****************************************************************************************
	//
	Array<T>& operator=(const Array<T>& ref) = default;

	// *****************************************************************************************
	//
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

	// *****************************************************************************************
	//
	void resize(int width, int height) {
		this->width = width;
		this->height = height;
		v.resize(width * height);
	}

	// *****************************************************************************************
	//
	void resize(int width, int height, const T& value) {
		this->width = width;
		this->height = height;
		v.resize(width * height, value);
	}

	// *****************************************************************************************
	//
	void clear() {
		v.clear();
	}

	// *****************************************************************************************
	//
	T* operator[](const int row) { return v.data() + row * width; }

	// *****************************************************************************************
	//
	const T* operator[](const int row) const { return v.data() + row * width; }

	// *****************************************************************************************
	//
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

// *****************************************************************************************
//
template <class T>
class LowerTriangularMatrix {
public:
	// *****************************************************************************************
	//
	size_t getSize() const { return size; }
	
	// *****************************************************************************************
	//
	LowerTriangularMatrix() : size(0) {}

	// *****************************************************************************************
	//
	LowerTriangularMatrix(size_t size) : size(size), data(size * (size - 1) / 2) {}

	// *****************************************************************************************
	//
	T* operator[](size_t i) { return data.data() + i * (i - 1) / 2; }

	// *****************************************************************************************
	//
	const T* operator[](size_t i) const { return data.data() + i * (i - 1) / 2; }

	// *****************************************************************************************
	//
	T* at(size_t i, size_t j) { return data[i * (i - 1) / 2  + j]; }

	// *****************************************************************************************
	//
	const T* at(size_t i, size_t j) const { return data[i * (i - 1) / 2 + j]; }

	// *****************************************************************************************
	//
	void resize(size_t size) {
		this->size = size;
		data.resize(size * (size - 1) / 2);
	}

	// *****************************************************************************************
	//
	void clear() {
		data.clear();
	}

	// *****************************************************************************************
	//
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

	void saveRow(size_t row_id, std::ofstream & file) {
		size_t offset = row_id * (row_id - 1) / 2;
		T * elem = data.data() + offset;

		for (size_t j = 0; j < row_id; ++j) {
			file << *elem++ << ',';
		}
	}

	void saveRowSparse(size_t row_id, std::ofstream & file) {
		size_t offset = row_id * (row_id - 1) / 2;
		T * elem = data.data() + offset;

		for (size_t j = 0; j < row_id; ++j, ++elem) {
			if (*elem > 0) {
				file << (j + 1) << ":" << *elem << ',';
			}
		}
	}

	int saveRow(size_t row, char* out) {
		size_t offset = row * (row - 1) / 2;
		return num2str(data.data() + offset, row, ',', out);
	}

	int saveRowSparse(size_t row_id, char* out) {
		size_t offset = row_id * (row_id - 1) / 2;
		T * elem = data.data() + offset;
		
		char *ptr = out;

		for (size_t j = 0; j < row_id; ++j, ++elem) {
			if (*elem > 0) {
				ptr += num2str(j + 1, ptr);
				*ptr++ = ':';
				ptr += num2str(*elem, ptr);
				*ptr++ = ',';
			}
		}
		
		return ptr - out;
	}


protected:
	size_t size;
	std::vector<T> data;

};


// *****************************************************************************************
//
template <typename T>
class heap {
public:
	heap() : currentEnd(data.begin()) {}
	heap(size_t maxSize) : data(maxSize), currentEnd(data.begin()) {}
	heap(const heap& ref) : data(ref.data) {
		currentEnd = data.begin();
	}

	void resize(size_t size) {
		data.resize(size);
		currentEnd = data.begin();
	}


	size_t size() const { return currentEnd - data.begin(); }

	bool push(const T& v) {
		bool stillBuilding = currentEnd != data.end();
		bool replaceMax = (currentEnd == data.end()) && (v < data.front());

		if (stillBuilding || replaceMax) {

			if (stillBuilding) {
				*currentEnd = v;
				++currentEnd;
			}

			if (replaceMax) {
				pop_heap(data.begin(), currentEnd);
				data.back() = v;
			}

			push_heap(data.begin(), currentEnd);
			return true;
		}

		return false;
	}

	std::vector<T>& getData() { return data;  }
	const std::vector<T>& getData() const { return data; }


protected:
	std::vector<T> data;
	typename std::vector<T>::iterator currentEnd;
};