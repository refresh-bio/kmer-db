#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cassert>

#include "conversion.h"
#include "../libs/refresh/sort/lib/pdqsort_par.h"

using namespace std;

// *****************************************************************************************
template<typename Item, typename Value, typename Score>
class Sampler
{
public:
	enum class strategy_t {none, best, random};

private:
	size_t no_samples;
	size_t max_items_per_sample;

	struct item_t
	{
		Item item;
		Value value;
		Score score;

		item_t(Item item, Value value, Score score) :
			item(item), value(value), score(score)
		{}

		item_t(const item_t&) = default;
		item_t(item_t&&) = default;
		item_t& operator=(const item_t&) = default;
		item_t& operator=(item_t&&) = default;
	};

	vector<vector<item_t>> data;
	vector<size_t> data_sizes;
	vector<mt19937_64> mts;
	strategy_t strategy{ strategy_t::none };

	static bool heap_comparer(const item_t& x, const item_t& y)
	{
		if (x.score != y.score)
			return x.score > y.score;
		return x.item < y.item;
	}

	void prepare_heap(size_t sample_id)
	{
		make_heap(data[sample_id].begin(), data[sample_id].end(), &(this->heap_comparer));
	}

	void select_best(size_t sample_id)
	{
		if (data[sample_id].back().score >= data[sample_id].front().score)			// if worse than min of heap, do not push into heap
		{
			push_heap(data[sample_id].begin(), data[sample_id].end(), &(this->heap_comparer));
			pop_heap(data[sample_id].begin(), data[sample_id].end(), &(this->heap_comparer));
		}
		data[sample_id].pop_back();
	}

	void select_random(size_t sample_id)
	{
		if (mts[sample_id]() % data_sizes[sample_id] == 0)			// remove the new one?
			;
		else
		{
			size_t id = mts[sample_id]() % max_items_per_sample;	// which to remove
			data[sample_id][id] = data[sample_id].back();
		}

		data[sample_id].pop_back();
	}

public:
	Sampler(size_t no_samples, size_t max_items_per_sample, strategy_t strategy) :
		no_samples(no_samples),
		max_items_per_sample(max_items_per_sample),
		strategy(strategy)
	{
		data.resize(no_samples);
		if (strategy == strategy_t::random)
		{
			data_sizes.resize(no_samples, 0);
			mts.resize(no_samples);
		}

		for (auto& x : data)
			x.reserve(max_items_per_sample + 1);
	}

	void add(size_t sample_id, Item item, Value value, Score score)
	{
		assert(sample_id < no_samples);

		data[sample_id].emplace_back(item, value, score);

		if(strategy == strategy_t::random)
			data_sizes[sample_id]++;

		if (strategy == strategy_t::best && data[sample_id].size() == max_items_per_sample)
			prepare_heap(sample_id);

		if (data[sample_id].size() <= max_items_per_sample)
			return;

		switch (strategy)
		{
		case strategy_t::best:
			select_best(sample_id);
			break;
		case strategy_t::random:
			select_random(sample_id);
			break;
		}
	}

	int saveRowSparse(size_t row_id, char* out, uint32_t idx_shift) 
	{
		auto out0 = out;

		refresh::sort::pdqsort(data[row_id].begin(), data[row_id].end(), [](const auto& x, const auto& y) {return x.item < y.item; });

		for (auto& x : data[row_id])
		{
			out += num2str(idx_shift + x.item + 1, out);
			*out++ = ':';
			out += num2str(x.value, out);
			*out++ = ',';
		}

		return (int) (out - out0);
	}
};