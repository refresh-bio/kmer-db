#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "conversion.h"
#include "../libs/refresh/pdqsort_par.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <unordered_map>
#include <atomic>
#include <future>

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

	const std::vector<T>& getData() const { return data; }
	std::vector<T>& getData() { return data; }
	
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

	int saveRow(size_t row_id, char* out) {
		size_t offset = row_id * (row_id - 1) / 2;
		return num2str(data.data() + offset, row_id, ',', out);
	}

	int saveRowSparse(size_t row_id, char* out) {
		size_t offset = row_id * (row_id - 1) / 2;
		return num2str_sparse(data.data() + offset, row_id, ',', out);
	}


protected:
	size_t size;
	std::vector<T> data;

};

#define SparseMatrix_CUSTOM_HASH_MAP

// *****************************************************************************************
//
template <class T>
class SparseMatrix {
#ifdef SparseMatrix_CUSTOM_HASH_MAP
	using hash_map_t = hash_map_lp<uint32_t, T>;
#else
	using hash_map_t = unordered_map<uint32_t, T>;
#endif
public:
	// *****************************************************************************************
	//
	size_t getSize() const { return size; }

	const std::vector<hash_map_t>& getData() const { return data; }
	std::vector<hash_map_t>& getData() { return data; }

	// *****************************************************************************************
	//
	SparseMatrix() : size(0) {}

	// *****************************************************************************************
	//
	SparseMatrix(size_t size) : size(size), data(size) {}

	// *****************************************************************************************
	//
	hash_map_t& operator[](size_t i) { return data[i]; }

	// *****************************************************************************************
	//
	const hash_map_t& operator[](size_t i) const { return data[i]; }

	// *****************************************************************************************
	//
	T* at(size_t i, size_t j) { return data[i][j]; }

	// *****************************************************************************************
	//
	const T* at(size_t i, size_t j) const { return data[i][j]; }

	// *****************************************************************************************
	//
	void resize(size_t size) {
		this->size = size;
		data.resize(size);
	}

	// *****************************************************************************************
	//
	void clear() {
		data.clear();
	}

	// *****************************************************************************************
	//
	void compact(uint32_t below, uint32_t above, uint32_t num_threads)
	{
		data_compacted.resize(data.size());

		if (tmp_rows.size() < num_threads)
			tmp_rows.resize(num_threads);

		atomic<uint32_t> row_id = 0;

		vector<future<void>> fut;
		fut.reserve(num_threads);

		for (uint32_t thread_id = 0; thread_id < num_threads; ++thread_id)
			fut.push_back(async([&, thread_id] {
				auto& my_row = tmp_rows[thread_id];
				while (true)
				{
					uint32_t i = row_id.fetch_add(1);

					if (i >= data.size())
						break;

					auto empty_value = data[i].empty_value;
					for (auto p = data[i].begin(); p != data[i].end(); ++p)
						if (p->val != empty_value && p->val < below && p->val > above)
							my_row.push_back(make_pair(p->key, p->val));

					refresh::sort::pdqsort_branchless(my_row.begin(), my_row.end());
					data_compacted[i].assign(my_row.begin(), my_row.end());
					my_row.clear();
				}

				my_row.shrink_to_fit();
			}
			));

		for (auto& f : fut)
			f.get();

/*		for (size_t i = 0; i < data.size(); ++i)
		{
			tmp_row.clear();

			auto empty_value = data[i].empty_value;
			for (auto p = data[i].begin(); p != data[i].end(); ++p)
				if (p->val != empty_value && p->val < below && p->val > above)
					tmp_row.push_back(make_pair(p->key, p->val));

			sort(tmp_row.begin(), tmp_row.end());
			data_compacted[i].assign(tmp_row.begin(), tmp_row.end());
		}

		tmp_row.clear();
		tmp_row.shrink_to_fit();*/

		data.clear();
		data.shrink_to_fit();
	}

	// *****************************************************************************************
	//
	void save(std::ofstream& file) {
		T* ptr = data.data();
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

/*	void saveRow(size_t row_id, std::ofstream& file) {
		size_t offset = row_id * (row_id - 1) / 2;
		T* elem = data.data() + offset;

		for (size_t j = 0; j < row_id; ++j) {
			file << *elem++ << ',';
		}
	}*/

	void saveRowSparse(size_t row_id, std::ofstream& file) {
#ifdef SparseMatrix_CUSTOM_HASH_MAP
		tmp_row.clear();

		auto empty_value = data[row_id].empty_value;
		for (auto p = data[row_id].begin(); p != data[row_id].end(); ++p)
			if (p->val != empty_value)
				tmp_row.push_back(make_pair(p->key, p->val));
#else
		tmp_row.assign(data[row_id].begin(), data[row_id].end());
#endif

//		sort(tmp_row.begin(), tmp_row.end());
		refresh::sort::pdqsort_branchless(tmp_row.begin(), tmp_row.end());

		for (const auto& x : tmp_row)
			if(x.first < row_id && x.second > 0)
				file << (x.first + 1) << ":" << x.second << ",";
	}

/*	int saveRow(size_t row_id, char* out) {
		size_t offset = row_id * (row_id - 1) / 2;
		return num2str(data.data() + offset, row_id, ',', out);
	}*/

	int saveRowSparse(size_t row_id, char* out) {
#ifdef SparseMatrix_CUSTOM_HASH_MAP
		tmp_row.clear();

		auto empty_value = data[row_id].empty_value;
		for (auto p = data[row_id].begin(); p != data[row_id].end(); ++p)
			if (p->val != empty_value)
				tmp_row.push_back(make_pair(p->key, p->val));
#else
		tmp_row.assign(data[row_id].begin(), data[row_id].end());
#endif

//		sort(tmp_row.begin(), tmp_row.end());
		refresh::sort::pdqsort_branchless(tmp_row.begin(), tmp_row.end());

		auto out0 = out;

		for (const auto& x : tmp_row)
			if (x.first < row_id && x.second > 0)
			{
				out += num2str(x.first + 1, out);
				*out++ = ':';
				out += num2str(x.second, out);
				*out++ = ',';
			}

		return out - out0;
	}

	int saveRowSparse(size_t row_id, char* out, uint32_t idx_shift) {
		auto out0 = out;

		for (const auto& x : data_compacted[row_id])
		{
			out += num2str(idx_shift + x.first + 1, out);
			*out++ = ':';
			out += num2str(x.second, out);
			*out++ = ',';
		}

		return out - out0;
	}


protected:
	size_t size;
	std::vector<hash_map_t> data;
	vector<vector<pair<uint32_t, T>>> data_compacted;
	vector<vector<pair<uint32_t, T>>> tmp_rows;
	vector<pair<uint32_t, T>> tmp_row;
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
