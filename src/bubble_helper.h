#pragma once

#include <vector>
#include <cinttypes>
#include <algorithm>

using namespace std;

class CBubbleHelper
{
public:
	struct bubble_adders_t
	{
		uint32_t num_kmers;
		vector<uint32_t>::iterator first;
		vector<uint32_t>::iterator last;

		bubble_adders_t(uint32_t num_kmers, vector<uint32_t>::iterator first, vector<uint32_t>::iterator last) :
			num_kmers(num_kmers), first(first), last(last)
		{}
	};

private:
	size_t bubble_size_thr;

	struct bubble_t
	{
		uint32_t num_kmers;
		vector<uint32_t> ids1;
		vector<uint32_t> ids2;

		template<typename Iter>
		bubble_t(uint32_t num_kmers, Iter first1, Iter last1, Iter first2, Iter last2) :
			num_kmers(num_kmers), ids1(first1, last1), ids2(first2, last2)
		{}

		template<typename Iter>
		bubble_t(uint32_t num_kmers, Iter first, Iter last) :
			num_kmers(num_kmers), ids1(first, last), ids2{}
		{}
	};

	vector<vector<uint32_t>> presence1;
	vector<vector<uint32_t>> presence2;
	
	vector<bubble_t> bubbles;

public:
	CBubbleHelper(size_t bubble_size_thr = 8000) :
		bubble_size_thr(bubble_size_thr)
	{}

	bool empty() const
	{
		return bubbles.empty();
	}

	void resize(size_t n)
	{
		clear();
		presence1.resize(n);
	}

	void resize(size_t n1, size_t n2)
	{
		clear();
		presence1.resize(n1);
		presence1.resize(n2);
	}

	void clear()
	{
		presence1.clear();
		presence2.clear();
		presence1.shrink_to_fit();
		presence2.shrink_to_fit();
	}

	bool is_bubble(size_t n)
	{
		return n >= bubble_size_thr;
	}

	bool is_bubble(size_t n1, size_t n2)
	{
		return n1 * n2 >= bubble_size_thr * bubble_size_thr;
	}

	template<typename Iter>
	void add(uint32_t num_kmers, Iter first, Iter last)
	{
		uint32_t bubble_id = (uint32_t) bubbles.size();

		uint32_t max_id = *(last - 1);
		if (presence1.size() <= max_id)
			presence1.resize(max_id + 1);

		for (auto p = first; p != last; ++p)
			presence1[*p].emplace_back(bubble_id);

		bubbles.emplace_back(num_kmers, first, last);
	}

	template<typename Iter>
	void add(uint32_t num_kmers, Iter first1, Iter last1, Iter first2, Iter last2)
	{
		uint32_t bubble_id = (uint32_t) bubbles.size();

		uint32_t max_id1 = *(last1 - 1);
		if (presence1.size() <= max_id1)
			presence1.resize(max_id1 + 1);

		// !!! TODO: Probably presence2 is unnecessary
		uint32_t max_id2 = *(last2 - 1);
		if (presence2.size() <= max_id2)
			presence2.resize(max_id2 + 1);

		for (auto p = first1; p != last1; ++p)
			presence1[*p].emplace_back(bubble_id);

		for (auto p = first2; p != last2; ++p)
			presence2[*p].emplace_back(bubble_id);

		bubbles.emplace_back(num_kmers, first1, last1, first2, last2);
	}

	bool get_adders(uint32_t id_query, vector<bubble_adders_t>& bubble_adders)
	{
		bubble_adders.clear();

		if (id_query >= presence1.size())
			return false;

		for (auto bid : presence1[id_query])
			if (bubbles[bid].ids2.empty())			// all2all-sp mode
			{
/*				auto p = upper_bound(bubbles[bid].ids1.begin(), bubbles[bid].ids1.end(), id_query);
				if (p != bubbles[bid].ids1.end())
					bubble_adders.emplace_back(bubbles[bid].num_kmers, p, bubbles[bid].ids1.end());*/

				auto p = lower_bound(bubbles[bid].ids1.begin(), bubbles[bid].ids1.end(), id_query);

				if (p != bubbles[bid].ids1.begin())
					bubble_adders.emplace_back(bubbles[bid].num_kmers, bubbles[bid].ids1.begin(), p);
			}
			else									// db2db-sp mode
			{
				bubble_adders.emplace_back(bubbles[bid].num_kmers, bubbles[bid].ids2.begin(), bubbles[bid].ids2.end());
			}

		return !bubble_adders.empty();
	}
};