#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "conversion.h"
#include "sampler.h"
#include "sparse_filters.h"
#include "../libs/refresh/sort/lib/pdqsort_par.h"
#include "bubble_helper.h"

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
	template <class Filter>
	void compact(const Filter& filter) {

		for (int i = 0; i < size; ++i) {
			T* row = (*this)[i];
			for (int j = 0; j < i; ++j) {
				if (!filter(row[j], i, j)) {
					row[j] = 0;
				}
			}
		}

	}

	// *****************************************************************************************
	//
	template <class Filter>
	void compact(const Filter& filter, uint32_t num_threads) {

		atomic<uint32_t> row_id = 0;

		vector<future<void>> fut;
		fut.reserve(num_threads);

		for (uint32_t thread_id = 0; thread_id < num_threads; ++thread_id) {
			fut.push_back(async([&, thread_id] {

				while (true) {
					uint32_t i = row_id.fetch_add(1);

					if (i >= data.size())
						break;

					T* row = (*this)[i];
					for (int j = 0; j < i; ++j) {
						if (!filter(row[j], i, j)) {
							row[j] = 0;
						}
					}
				}
			}));
		}

		for (auto& f : fut) {
			f.get();
		}
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
	void compact(uint32_t mini, uint32_t maxi, uint32_t num_threads)
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
						if (p->val != empty_value && p->val >= mini && p->val <= maxi)
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
	template <class Filter>
	void compact2(const Filter& filter, uint32_t num_threads, CBubbleHelper &bubbles)
	{
		data_compacted.resize(data.size());

		if (tmp_rows.size() < num_threads)
			tmp_rows.resize(num_threads);

		atomic<uint32_t> row_id = 0;

		vector<future<void>> fut;
		fut.reserve(num_threads);

		for (uint32_t thread_id = 0; thread_id < num_threads; ++thread_id)
			fut.push_back(async([&, thread_id] {
			vector<CBubbleHelper::bubble_adders_t> bubble_adders;

			auto& my_row = tmp_rows[thread_id];
			while (true)
			{
				uint32_t i = row_id.fetch_add(1);

				if (i >= data.size())
					break;

				if (!bubbles.empty())
				{
					bubbles.get_adders(i, bubble_adders);
					for (auto adder : bubble_adders)
						for (auto p = adder.first; p != adder.last; ++p)
							data[i][*p] += adder.num_kmers;
				}

				auto empty_value = data[i].empty_value;
				for (auto p = data[i].begin(); p != data[i].end(); ++p)
					if (p->val != empty_value && filter(p->val, i, p->key))
						my_row.push_back(make_pair(p->key, p->val));

				data[i].release();

				refresh::sort::pdqsort_branchless(my_row.begin(), my_row.end());
				data_compacted[i].assign(my_row.begin(), my_row.end());
				my_row.clear();
			}

			my_row.shrink_to_fit();
				}
			));

		for (auto& f : fut)
			f.get();


		data.clear();
		data.shrink_to_fit();
	}


	// *****************************************************************************************
	void add_to_sampler(
		CombinedFilter<uint32_t> &filter, 
		Sampler<uint32_t, uint32_t, double> &sampler, 
		const metric_fun_t &sampling_criterion, 
		const std::vector<uint32_t>& num_row_kmers, 
		const std::vector<uint32_t>& num_col_kmers, 
		const uint32_t idx_row_shift,
		const uint32_t idx_col_shift,
		const uint32_t k_len,
		const uint32_t num_threads,
		CBubbleHelper &bubbles)
	{
		atomic<uint32_t> row_id = 0;

		vector<uint32_t> max_col_ids(num_threads, 0);

		vector<future<void>> fut;
		fut.reserve(num_threads);

		vector<vector<pair<uint32_t, T>>> accepted_items(data.size());

		for (uint32_t thread_id = 0; thread_id < num_threads; ++thread_id)
			fut.push_back(async([&, thread_id] {
			auto local_sampling_criterion = sampling_criterion ? sampling_criterion : [](num_kmers_t common, num_kmers_t seq1, num_kmers_t seq2, int k_len) ->double {return 1; };
			uint32_t my_max_col_id = 0;
			vector<CBubbleHelper::bubble_adders_t> bubble_adders;

			while (true)
			{
				uint32_t i = row_id.fetch_add(1);

				if (i >= data.size())
					break;

				if (!bubbles.empty())
				{
					bubbles.get_adders(i, bubble_adders);
					for (auto adder : bubble_adders)
						for (auto p = adder.first; p != adder.last; ++p)
							data[i][*p] += adder.num_kmers;
				}

				auto empty_value = data[i].empty_value;
				for (auto p = data[i].begin(); p != data[i].end(); ++p)
					if (p->val != empty_value && filter(p->val, i, p->key))
					{
						sampler.add(i + idx_row_shift, p->key + idx_col_shift, p->val, (double)local_sampling_criterion(p->val, num_row_kmers[i], num_col_kmers[p->key], k_len));
						accepted_items[i].emplace_back(p->key, p->val);

						if (p->key > my_max_col_id)
							my_max_col_id = p->key;
					}

				data[i].release();
			}

			max_col_ids[thread_id] = my_max_col_id;
				}
			));

		for (auto& f : fut)
			f.get();
		fut.clear();

		// Transpose
		uint32_t max_col_id = *max_element(max_col_ids.begin(), max_col_ids.end());
		vector<vector<pair<uint32_t, T>>> data_transposed(max_col_id + 1);

		for (uint32_t i = 0; i < accepted_items.size(); ++i)
			for (auto& x : accepted_items[i])
				data_transposed[x.first].emplace_back(i, x.second);

		atomic<uint32_t> tr_row_id = 0;

		for (uint32_t thread_id = 0; thread_id < num_threads; ++thread_id)
			fut.push_back(async([&, thread_id] {
			auto local_sampling_criterion = sampling_criterion ? sampling_criterion : [](num_kmers_t common, num_kmers_t seq1, num_kmers_t seq2, int k_len) ->double {return 1; };

			while (true)
			{
				uint32_t i = tr_row_id.fetch_add(1);

				if (i >= data_transposed.size())
					break;

				for (auto &x : data_transposed[i])
					sampler.add(i + idx_col_shift, x.first + idx_row_shift, x.second, (double)local_sampling_criterion(x.second, num_row_kmers[x.first], num_col_kmers[i], k_len));
			}
				}
			));

		for (auto& f : fut)
			f.get();
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

		return (int) (out - out0);
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

		return (int) (out - out0);
	}

	size_t getNoInRow(size_t row_id)	const
	{
		return data_compacted[row_id].size();
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
