#pragma once

#include "record.h"
#include "exceptions.h"
#include "small_sort.h"
#include "defs.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <mutex>
#include <thread>
#include <cstring>
#include <condition_variable>
#include <memory>
#include <emmintrin.h>
#include <immintrin.h>

namespace raduls
{
	//config
	const uint32_t MAX_REC_SIZE_IN_BYTES = 32;

	constexpr uint32 ALIGNMENT = 0x100;
	constexpr int32 BUFFER_WIDTHS[] = { -1, 32, 16, 16, 8, 8, 4, 8, 4 };
	const uint64 small_sort_thresholds[] = { 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384 };	
	const uint64 wide_small_sort_thresholds[] = { 64, 48, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 };

	const uint32 FIRST_PASS_THREADS_MULTIPLIER = 8;
	const uint64 SMALL_RADIX_THRESHOLD = 256;

	static_assert(MAX_REC_SIZE_IN_BYTES / 8 <= (sizeof(BUFFER_WIDTHS) / sizeof(BUFFER_WIDTHS[0]) - 1), "BUFFER_WIDTHS must be extended");
	static_assert(MAX_REC_SIZE_IN_BYTES % 8 == 0, "MAX_REC SIZE MUST BE A MULTIPLE OF 8");

#ifdef __GNUG__
#if __GNUC__ < 6
	static_assert(false, "gcc 6 or higher is required.");
#endif
#endif
#ifdef _MSC_VER
#if _MSC_VER < 1900
	static_assert(false, "Visual Studio 2015 or higher is required.");
#endif
#endif

	inline bool check_narrowing(uint64 x, uint64 y) { return x < 16 * y; }
	class CRangeQueue
	{
		std::vector<std::tuple<uint64, uint64, uint32>> range_queue;
		std::mutex mtx;
		uint32 cur_idx;
		bool done;
	public:
		CRangeQueue(uint32 parts, uint64 num_rec)
		{
			uint64 delta = num_rec / parts;
			uint64 N1 = 0, N2 = 0;
			uint64 smallest_fraction = 8;
			uint64 start_size = delta / smallest_fraction;
			uint64 step = (2 * smallest_fraction - 2) * delta / smallest_fraction / parts;
			uint64 cur_delta = start_size;

			for (uint32 i = 0; i < parts; ++i)
			{
				N2 = N1 + cur_delta;
				cur_delta += step;
				if (i == parts - 1)
					N2 = num_rec;
				range_queue.emplace_back(N1, N2, parts - i - 1);
				N1 = N2;
			}
			std::reverse(range_queue.begin(), range_queue.end());

			cur_idx = 0;
			if (parts)
				done = false;
		}

		bool get(uint64 &n1, uint64 &n2, uint32 &part_id)
		{
			std::lock_guard<std::mutex> lck(mtx);

			if (!done)
			{
				std::tie(n1, n2, part_id) = range_queue[cur_idx++];
				if (cur_idx == range_queue.size())
					done = true;
				return true;
			}
			return false;
		}

		void reset_indices()
		{
			cur_idx = 0;
			done = false;
		}
	};

	template<typename RECORD_T>
	struct CRadixMSDTaskskDesc
	{
		RECORD_T *data, *tmp;
		uint64 n_recs;
		uint32 byte;
		bool is_narrow;

		CRadixMSDTaskskDesc(RECORD_T* data, RECORD_T *tmp, uint64 n_recs, uint32 byte, bool is_narrow) :
			data(data),
			tmp(tmp),
			n_recs(n_recs),
			byte(byte),
			is_narrow(is_narrow)
		{ }

		bool operator<(const CRadixMSDTaskskDesc<RECORD_T> &rhs) const
		{
			return this->n_recs < rhs.n_recs;
		}
	};

	template<typename RECORD_T>
	class CRadixMSDTaskQueue
	{
		std::priority_queue<CRadixMSDTaskskDesc<RECORD_T>> tasks;
		std::condition_variable cv_pop;
		std::mutex mtx;
		uint64 tasks_in_progress = 0;

	public:
		void push(RECORD_T* data, RECORD_T* tmp, uint64 n_recs, uint32 byte, bool is_narrow)
		{
			std::lock_guard<std::mutex> lck(mtx);
			tasks_in_progress++;
			tasks.emplace(data, tmp, n_recs, byte, is_narrow);
			if (tasks.size() == 1) //was empty
				cv_pop.notify_all();
		}

		bool pop(RECORD_T* &data, RECORD_T* &tmp, uint64 &n_recs, uint32 &byte, bool &is_narrow)
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this] {return tasks.size() || !tasks_in_progress; });
			if (!tasks_in_progress)
				return false;

			data = tasks.top().data;
			tmp = tasks.top().tmp;
			n_recs = tasks.top().n_recs;
			byte = tasks.top().byte;
			is_narrow = tasks.top().is_narrow;
			tasks.pop();
			return true;
		}

		void notify_task_finished()
		{
			std::lock_guard<std::mutex> lck(mtx);
			--tasks_in_progress;
			if (!tasks_in_progress)
				cv_pop.notify_all();
		}
	};

	// 64b copy function
	// size - in 8B words (determined during execution)
	// dest and src must be aligned to 8B
	inline void IntrCopy64fun(void *_dest, void *_src, uint32_t size)
	{
		__int64* dest = (__int64 *)_dest;
		__int64* src = (__int64 *)_src;

		for (unsigned i = 0; i < size; ++i)			
			_mm_stream_si64(dest + i, src[i]);	
	}

	// 64bit copy function
	// SIZE - in 8B words
	template <unsigned SIZE> struct IntrCopy64
	{
		static inline void Copy(void *_dest, void *_src)
		{
			__int64* dest = (__int64*)_dest;
			__int64* src = (__int64*)_src;

			for (unsigned i = 0; i < SIZE; ++i)
				_mm_stream_si64(dest + i, src[i]);
		}
	};


	template <unsigned SIZE, unsigned MODE> struct IntrCopy128
	{
	
	};

	// 128bit copy function
	// SIZE - in 16B words
	// dest - aligned to 16B
	// src  - aligned to 16B
	template <unsigned SIZE> struct IntrCopy128<SIZE, 1>
	{
		static inline void Copy(void *_dest, void *_src)
		{
			__m128i *dest = (__m128i *) _dest;
			__m128i *src = (__m128i *) _src;

			for (unsigned i = 0; i < SIZE; ++i)
				_mm_stream_si128(dest + i, _mm_load_si128(src + i));
		}
	};


	// 128bit copy function
	// SIZE - in 16B words
	// dest - aligned to 8B
	// src  - aligned to 16B
	template <unsigned SIZE> struct IntrCopy128<SIZE, 0>
	{
		static inline void Copy(void *dest, void *src)
		{
			if ((uint64_t)dest % 16)	// if only 8B aligned use 64b copy
				IntrCopy64<SIZE * 2>::Copy(dest, src);
			else // if 16B aligned use 128b copy
				IntrCopy128<SIZE, 1>::Copy(dest, src);
		}
	};


	template<typename RECORD_T, typename COUNTER_TYPE>
	FORCE_INLINE void BuildHisto(COUNTER_TYPE* histo, uint64 n, uint8_t* &ptr)
	{
		switch (n % 4)
		{
		case 3:
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		case 2:
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		case 1:
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		}

		auto n_iters = n / 4;
		for (uint64 i = 0; i < n_iters; ++i)
		{
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		}
	}

	template<typename RECORD_T, typename COUNTER_TYPE>
	FORCE_INLINE void SimpleScatter(RECORD_T* src, RECORD_T* tmp, COUNTER_TYPE* histo, uint64 n, uint8_t* &ptr)
	{
		switch (n % 4)
		{
		case 3:
			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		case 2:
			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		case 1:
			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		}

		for (uint64 i = n % 4; i < n; i += 4)
		{
			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);

			tmp[histo[*ptr]] = *src++;
			histo[*ptr]++;
			ptr += sizeof(RECORD_T);
		}
	}

	template<typename RECORD_T, typename COUNTER_TYPE,
		uint32 BUFFER_WIDTH, uint32 BUFFER_WIDTH_IN_128BIT_WORDS>
		FORCE_INLINE void BufferedScatterStep(RECORD_T* &src, RECORD_T* tmp,
			COUNTER_TYPE* histo, COUNTER_TYPE* copy_histo, uint8_t* &ptr, uint8_t &byteValue,
			RECORD_T* buffer, int &index_x, uint8_t* first_store)
	{
		byteValue = *ptr;
		index_x = histo[byteValue] % BUFFER_WIDTH;
		buffer[byteValue * BUFFER_WIDTH + index_x] = *src++;
		histo[byteValue]++;
		if (index_x == (BUFFER_WIDTH - 1))
		{
			if (first_store[byteValue])
			{
				first_store[byteValue] = false;
				int64 offset = copy_histo[byteValue] % BUFFER_WIDTH;
				IntrCopy64fun(&tmp[histo[byteValue] - BUFFER_WIDTH + offset], &buffer[byteValue * BUFFER_WIDTH + offset], RECORD_T::RECORD_SIZE * (BUFFER_WIDTH - offset));
			}
			else
				IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, 1>::Copy(&tmp[histo[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
		}
		ptr += sizeof(RECORD_T);
	}


	template<typename RECORD_T, typename COUNTER_TYPE, uint32 BUFFER_WIDTH>
	FORCE_INLINE void BufferedScattterCorrectionStep(RECORD_T* tmp, COUNTER_TYPE* histo, COUNTER_TYPE* copy_histo, RECORD_T* buffer)
	{
		int64 elemInBuffer, index_stop, index_start, elemWrittenIntoBuffer;
		for (uint32 i = 0; i < 256; ++i)
		{
			index_stop = histo[i] % BUFFER_WIDTH;
			index_start = copy_histo[i] % BUFFER_WIDTH;
			elemWrittenIntoBuffer = histo[i] - copy_histo[i];

			if ((index_stop - elemWrittenIntoBuffer) <= 0)
				elemInBuffer = index_stop;
			else
				elemInBuffer = index_stop - index_start;

			if (elemInBuffer != 0)
				IntrCopy64fun(&tmp[histo[i] - elemInBuffer],
					&buffer[i * BUFFER_WIDTH + (histo[i] - elemInBuffer) % BUFFER_WIDTH], elemInBuffer * sizeof(RECORD_T) / 8);
		}
	}

	template<typename RECORD_T, typename COUNTER_TYPE>
	void FirstPassStage1(RECORD_T* data, std::vector<COUNTER_TYPE[256]>& histos, uint32 byte, CRangeQueue& rq)
	{
		alignas(ALIGNMENT) COUNTER_TYPE myHisto[256] = {}; 
		uint64 idx1, idx2;
		uint32 part_id;

		while (rq.get(idx1, idx2, part_id))
		{
			std::memset(myHisto, 0, sizeof(myHisto));

			auto ptr = reinterpret_cast<uint8_t*>(data + idx1) + byte;
			uint64 n = idx2 - idx1;

			BuildHisto<RECORD_T, COUNTER_TYPE>(myHisto, n, ptr);

			for (uint32 i = 0; i < 256; ++i)
				histos[part_id][i] = myHisto[i];
		}
	}
	template<typename RECORD_T, typename COUNTER_TYPE>
	void BigBinsScatter(RECORD_T* data, RECORD_T* tmp,
		uint32 byte, std::vector<COUNTER_TYPE[256]>& histos,
		std::vector<uchar*>& buffers, std::vector<COUNTER_TYPE[256]> &threads_histos,
		CRangeQueue& rq)
	{
		alignas(ALIGNMENT)COUNTER_TYPE myHisto[256];

		uint64 idx1, idx2;
		uint32 part_id;

		uint8_t* ptr;
		uint64 n;
		RECORD_T* src;

		constexpr uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(RECORD_T) / 8];
		constexpr uint32 BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(RECORD_T) / 16;		

		uint8_t byteValue = 0;
		int index_x = 0;

		while (rq.get(idx1, idx2, part_id))
		{
			for (int i = 0; i < 256; ++i)
				myHisto[i] = histos[part_id][i];

			auto copy_histo = histos[part_id];

			ptr = reinterpret_cast<uint8_t*>(data + idx1) + byte;
			n = idx2 - idx1;

			auto buffer = reinterpret_cast<RECORD_T*>(buffers[part_id]);

			src = data + idx1;
			byteValue = 0;
			index_x = 0;

			uint8_t first_store[256];
			std::fill_n(first_store, 256, true);

			switch (n % 4)
			{
			case 3:
				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);
			case 2:
				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);
			case 1:
				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);
			}

			for (uint64 i = n % 4; i < n; i += 4)
			{
				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);

				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);

				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);

				BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
					(src, tmp, myHisto, copy_histo, ptr, byteValue, buffer, index_x, first_store);
			}

			for (uint32 i = 0; i < 256; ++i)
				threads_histos[part_id][i] = myHisto[i];
		}
	}

	template<typename RECORD_T, typename COUNTER_TYPE>
	void FirstPassStage2(RECORD_T* data, RECORD_T* tmp, 
		uint32 byte, std::vector<COUNTER_TYPE[256]>& histos,
		std::vector<uchar*>& buffers, std::vector<COUNTER_TYPE[256]> &threads_histos,
		CRangeQueue& rq)
	{
		alignas(ALIGNMENT) COUNTER_TYPE myHisto[256];

		uint64 idx1, idx2;
		uint32 part_id;

		uint8_t* ptr;
		uint64 n;
		RECORD_T* src;

		constexpr uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(RECORD_T) / 8];
		constexpr uint32 BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(RECORD_T) / 16;
		constexpr uint32 BUFFER_16B_ALIGNED = 1; 

		uint8_t byteValue = 0;
		int index_x = 0;

		while (rq.get(idx1, idx2, part_id))
		{
			for (int i = 0; i < 256; ++i)
				myHisto[i] = histos[part_id][i];

			ptr = reinterpret_cast<uint8_t*>(data + idx1) + byte;
			n = idx2 - idx1;

			auto buffer = reinterpret_cast<RECORD_T*>(buffers[part_id]);

			src = data + idx1;
			byteValue = 0;
			index_x = 0;

			switch (n % 4)
			{
			case 3:
				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 3];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);
			case 2:
				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 2];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);
			case 1:
				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[(n % 4) - 1];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);
			}

			for (uint64 i = n % 4; i < n; i += 4)
			{
				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[i];
				myHisto[byteValue]++;
				if(index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);

				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 1];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);

				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 2];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);

				byteValue = *ptr;
				index_x = myHisto[byteValue] % BUFFER_WIDTH;
				buffer[byteValue * BUFFER_WIDTH + index_x] = src[i + 3];
				myHisto[byteValue]++;
				if (index_x == (BUFFER_WIDTH - 1))
					IntrCopy128<BUFFER_WIDTH_IN_128BIT_WORDS, BUFFER_16B_ALIGNED>::Copy(&tmp[myHisto[byteValue] - BUFFER_WIDTH], &buffer[byteValue * BUFFER_WIDTH]);
				ptr += sizeof(RECORD_T);
			}

			for (uint32 i = 0; i < 256; ++i)
				threads_histos[part_id][i] = myHisto[i];
		}
	}

	template<typename RECORD_T, typename COUNTER_TYPE>
	void FirstPassStage3(RECORD_T* tmp, 
		std::vector<COUNTER_TYPE[256]>& histos,
		std::vector<uchar*>& buffers, std::vector<COUNTER_TYPE[256]> &threads_histos,
		CRangeQueue& rq)
	{
		uint64 idx1, idx2;
		uint32 part_id;
		alignas(ALIGNMENT)COUNTER_TYPE myHisto[256];
		
		const uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(RECORD_T) / 8];

		while (rq.get(idx1, idx2, part_id))
		{
			auto buffer = reinterpret_cast<RECORD_T*>(buffers[part_id]);

			for (int i = 0; i < 256; ++i)
				myHisto[i] = threads_histos[part_id][i];

			BufferedScattterCorrectionStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH>(tmp, myHisto, histos[part_id], buffer);
		}
	}

	template<typename RECORD_T>
	void SmallSortDispatch(RECORD_T* data, RECORD_T* tmp, uint64 size)
	{
		small_sort::HybridSmallSort<RECORD_T>{}(data, size);
		//std::sort(data, data + size);		
	}

	template<typename RECORD_T, typename COUNTER_TYPE>
	class CRadixSorterMSD
	{
		CRadixMSDTaskQueue<RECORD_T>& tasks_queue;
		uint64 use_queue_min_recs = 0;
		uchar* _buffer;
		void Sort(RECORD_T* data, RECORD_T* tmp, uint64 n_recs, uint32 byte, bool is_narrow)
		{
			auto ptr = reinterpret_cast<uint8_t*>(data) + byte;
			alignas(ALIGNMENT) COUNTER_TYPE globalHisto[256] = {};
			alignas(ALIGNMENT) COUNTER_TYPE copy_globalHisto[257];
			uint64 largest_bin_size = 0;

			BuildHisto<RECORD_T, COUNTER_TYPE>(globalHisto, n_recs, ptr);

			COUNTER_TYPE prevSum = 0;
			for (int i = 0; i < 256; ++i)
			{
				COUNTER_TYPE tmp = globalHisto[i];
				globalHisto[i] = prevSum;
				copy_globalHisto[i] = prevSum;
				prevSum += tmp;

				if (tmp > largest_bin_size)
					largest_bin_size = tmp;
			}
			copy_globalHisto[256] = static_cast<COUNTER_TYPE>(n_recs);

			auto src = data;
			ptr = reinterpret_cast<uint8_t*>(data) + byte;

			//if fits in L2 cache
			if (n_recs * sizeof(RECORD_T) < (1ull << 16))			
				SimpleScatter(src, tmp, globalHisto, n_recs, ptr);			
			else
			{
				constexpr uint32 BUFFER_WIDTH = BUFFER_WIDTHS[sizeof(RECORD_T) / 8];
				constexpr uint32 BUFFER_WIDTH_IN_128BIT_WORDS = BUFFER_WIDTH * sizeof(RECORD_T) / 16;

				auto buffer = reinterpret_cast<RECORD_T*>(_buffer);
			
				uint8_t byteValue = 0;
				int index_x = 0;

				uint8_t first_store[256];
				std::fill_n(first_store, 256, true);
			
				// Move back tmp pointer - to be aligned to 64B
				uint64 tmp_moved_by = 0;
				while ((uint64)tmp % 64)
					tmp--, tmp_moved_by++;
				
				// Update histograms - they must point correct places after tmp alignment
				for (int i = 0; i < 256; ++i)
				{
					globalHisto[i] += tmp_moved_by;
					copy_globalHisto[i] += tmp_moved_by;
 				}
				copy_globalHisto[256] += tmp_moved_by;

				switch (n_recs % 4)
				{
				case 3:
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
				case 2:
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
				case 1:
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
				}
				
				for (uint64 i = n_recs % 4; i < n_recs; i += 4)
				{
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
					BufferedScatterStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH, BUFFER_WIDTH_IN_128BIT_WORDS>
						(src, tmp, globalHisto, copy_globalHisto, ptr, byteValue, buffer, index_x, first_store);
				}

				BufferedScattterCorrectionStep<RECORD_T, COUNTER_TYPE, BUFFER_WIDTH>(tmp, globalHisto, copy_globalHisto, buffer);

				// Bring back tmp pointer to the original (unaligned) value
				tmp += tmp_moved_by;
				for (int i = 0; i < 256; ++i)
				{
					globalHisto[i] -= tmp_moved_by;
					copy_globalHisto[i] -= tmp_moved_by;
				}
				copy_globalHisto[256] -= tmp_moved_by;
			}

			if (byte > 0)
			{
				bool must_copy_tmp = byte % 2;

				uint64 narrow_small_sort_threshold = small_sort_thresholds[RECORD_T::RECORD_SIZE];
				uint64 wide_small_sort_threshold = wide_small_sort_thresholds[RECORD_T::RECORD_SIZE];

				if (largest_bin_size <= wide_small_sort_threshold || (is_narrow && largest_bin_size <= narrow_small_sort_threshold))
				{
					// All bins are small
					for (int i = 0; i < 256; ++i)
					{
						uint64 new_n = copy_globalHisto[i + 1] - copy_globalHisto[i];
						if (new_n > 1)
							SmallSortDispatch<RECORD_T>(tmp + copy_globalHisto[i], data + copy_globalHisto[i], new_n);
					}

					IntrCopy64fun(data, tmp, n_recs * sizeof(RECORD_T) / 8);					
				}
				else
					for (int i = 0; i < 256; ++i)
					{
						uint64 new_n = copy_globalHisto[i + 1] - copy_globalHisto[i];

						if (new_n <= wide_small_sort_threshold || (is_narrow && new_n <= narrow_small_sort_threshold))
						{
							if (new_n > 1)
								SmallSortDispatch(tmp + copy_globalHisto[i], data + copy_globalHisto[i], new_n);
							if (must_copy_tmp && new_n)
								for (COUNTER_TYPE j = copy_globalHisto[i]; j < copy_globalHisto[i] + static_cast<COUNTER_TYPE>(new_n); ++j)
									data[j] = tmp[j];
						}
						else
						{							
							if (new_n >= use_queue_min_recs)
								tasks_queue.push(tmp + copy_globalHisto[i], data + copy_globalHisto[i], new_n, byte - 1, is_narrow);
							else
								if (new_n < SMALL_RADIX_THRESHOLD)
									SmallRadixSort(tmp + copy_globalHisto[i], data + copy_globalHisto[i], new_n, byte - 1);
								else
									Sort(tmp + copy_globalHisto[i], data + copy_globalHisto[i], new_n, byte - 1, check_narrowing(n_recs, new_n));
						}
					}
			}			
			_mm_sfence();
		}

		void SmallRadixSort(RECORD_T* data, RECORD_T* tmp, uint64 n_recs, uint32 byte)
		{			
			auto ptr = reinterpret_cast<uint8_t*>(data) + byte;
			alignas(ALIGNMENT) uint32 globalHisto[256] = {};
			alignas(ALIGNMENT) uint32 copy_globalHisto[257];
			int to_sort[256];		// bins to sort (at least 2 elements in bin)
			int idx_to_sort = 0;
			bool must_copy_tmp = byte % 2;

			BuildHisto<RECORD_T, uint32>(globalHisto, n_recs, ptr);

			uint32 prevSum = 0;
			for (int i = 0; i < 256; ++i)
			{
				uint32 n_elems = globalHisto[i];
				globalHisto[i] = prevSum;
				copy_globalHisto[i] = prevSum;
				prevSum += n_elems;

				to_sort[idx_to_sort] = i;
				idx_to_sort += n_elems > 1;
			}
			copy_globalHisto[256] = static_cast<COUNTER_TYPE>(n_recs);

			auto src = data;
			ptr = reinterpret_cast<uint8_t*>(data) + byte;

			SimpleScatter(src, tmp, globalHisto, n_recs, ptr);

			if (byte > 0)
			{
				for (int ii = 0; ii < idx_to_sort; ++ii)
				{
					int i = to_sort[ii];
					uint64 new_n = copy_globalHisto[i + 1] - copy_globalHisto[i];

					SmallSortDispatch(tmp + copy_globalHisto[i], tmp + copy_globalHisto[i], new_n);
				}
				if (must_copy_tmp)
					IntrCopy64fun(data, tmp, n_recs * sizeof(RECORD_T) / 8);
			}
			_mm_sfence();
		}
	public:
		CRadixSorterMSD(CRadixMSDTaskQueue<RECORD_T>& tasks_queue, uint64 use_queue_min_recs, uchar* _buffer)
			:
			tasks_queue(tasks_queue),
			use_queue_min_recs(use_queue_min_recs),
			_buffer(_buffer)
		{}

		void operator()()
		{
			RECORD_T *data, *tmp;
			uint64 n_recs;
			uint32 byte;
			bool is_narrow;

			while (tasks_queue.pop(data, tmp, n_recs, byte, is_narrow))
			{
				Sort(data, tmp, n_recs, byte, is_narrow);
				tasks_queue.notify_task_finished();
			}
		}
	};

	template<typename RECORD_T>
	class CRaduls
	{
		uint32 n_threads;			

	public:
		CRaduls(uint32 n_threads):
			n_threads(n_threads)
		{			
		}

		template<typename COUNTER_TYPE>
		void Sort(RECORD_T* data, RECORD_T* tmp, uint64 n_recs, uint32 byte,
			uint32 n_threads, bool is_first_level, uint64 is_big_threshold, uint64 n_total_recs)
		{		
			uint64 current_small_sort_threshold = small_sort_thresholds[RECORD_T::RECORD_SIZE];

			if (n_recs <= current_small_sort_threshold)
			{
				SmallSortDispatch(data, tmp, n_recs);
				if (byte % 2 == 0)
					for (uint64 j = 0; j < n_recs; ++j)
						tmp[j] = data[j];
				return;
			}

			//stage 1
			const auto n_parts = FIRST_PASS_THREADS_MULTIPLIER * n_threads;
			CRangeQueue range_queue(n_parts, n_recs);
			std::vector<std::thread> threads;
			std::vector<COUNTER_TYPE[256]> histos(n_parts); 
			alignas(ALIGNMENT) COUNTER_TYPE globalHisto[257] = {};

			
			for (uint32 th_id = 0; th_id < n_threads; ++th_id)
				threads.emplace_back(FirstPassStage1<RECORD_T, COUNTER_TYPE>,
					data, std::ref(histos), byte, std::ref(range_queue));			

			for (auto& th : threads)
				th.join();
			threads.clear();


			// ***** collecting counters
			for (int i = 0; i < 256; ++i)
			{
				COUNTER_TYPE prevSum = 0;
				for (uint32 n = 0; n < n_parts; ++n)
				{
					auto tmp = histos[n][i];
					histos[n][i] = prevSum;
					prevSum += tmp;
				}
				globalHisto[i] = prevSum;
			}

			COUNTER_TYPE prevSum = 0;
			for (int i = 0; i < 256; ++i)
			{
				COUNTER_TYPE tmp = globalHisto[i];
				globalHisto[i] = prevSum;
				prevSum += tmp;
			}

			for (uint32 n = 0; n < n_parts; ++n)
				for (int i = 0; i < 256; ++i)
					histos[n][i] += globalHisto[i];

			//stage 2
			range_queue.reset_indices();
			
			uint64 single_part_size = 256 * BUFFER_WIDTHS[sizeof(RECORD_T)/8] * sizeof(RECORD_T);
			auto _raw_buffers = std::make_unique<uchar[]>(single_part_size * n_parts + ALIGNMENT);
			auto s = _raw_buffers.get();
			while ((uint64)s % ALIGNMENT)
				++s;
			std::vector<uchar*> buffers(n_parts);
			for (uint32 i = 0; i < n_parts; ++i)
				buffers[i] = s + single_part_size * i;
			
			std::vector<COUNTER_TYPE[256]> threads_histos(n_parts); 

			uint64 tmp_moved_by = 0;
			if (!is_first_level)
			{
				while ((uint64)tmp % 64)
					tmp--, tmp_moved_by++;

				for (uint32 n = 0; n < n_parts; ++n)
					for (int i = 0; i < 256; ++i)
						histos[n][i] += tmp_moved_by;
			}

			auto fun = is_first_level ? FirstPassStage2<RECORD_T, COUNTER_TYPE> : BigBinsScatter<RECORD_T, COUNTER_TYPE>;

			for (uint32 th_id = 0; th_id < n_threads; ++th_id)
				threads.emplace_back(fun,
					data, tmp, byte,
					std::ref(histos), std::ref(buffers), std::ref(threads_histos),
					std::ref(range_queue));

			for (auto& th : threads)
				th.join();
			threads.clear();

			//stage 3
			range_queue.reset_indices();
			for (uint32 th_id = 0; th_id < n_threads; ++th_id)
				threads.emplace_back(FirstPassStage3<RECORD_T, COUNTER_TYPE>,
					tmp, std::ref(histos), std::ref(buffers), std::ref(threads_histos),
					std::ref(range_queue));

			for (auto& th : threads)
				th.join();
			threads.clear();

			if (!is_first_level)
			{
				tmp += tmp_moved_by;
				for (uint32 n = 0; n < n_parts; ++n)
					for (int i = 0; i < 256; ++i)
						histos[n][i] -= tmp_moved_by;
			}

			if (byte > 0)
			{
				CRadixMSDTaskQueue<RECORD_T> tasks_queue;

				auto data_ptr = data;
				auto ptr = tmp;

				std::vector<std::tuple<RECORD_T*, RECORD_T*, uint64>> big_bins;
				uint64 n_recs_in_big_bins = 0;

				globalHisto[256] = n_recs;
				for (uint32 i = 1; i < 257; ++i)
				{
					auto n = static_cast<uint64>(globalHisto[i] - globalHisto[i - 1]);
					if (n > 0)
					{
						if (n > is_big_threshold)
						{
							if(!is_first_level)
								Sort<COUNTER_TYPE>(ptr, data_ptr, n, byte - 1, n_threads, false, is_big_threshold, n_total_recs);
							else
								big_bins.emplace_back(ptr, data_ptr, n),
								n_recs_in_big_bins += n;
						}
						else
							tasks_queue.push(ptr, data_ptr, n, byte - 1, check_narrowing(n_recs, n));
					}

					ptr += n;
					data_ptr += n;
				}

				std::sort(big_bins.begin(), big_bins.end(), [](const auto& x, const auto& y) {return std::get<2>(x) > std::get<2>(y); });

				auto n_threads_for_big_bins = std::min(n_threads, static_cast<uint32>(ceil(n_threads * n_recs_in_big_bins * 5.0 / (4 * n_total_recs))));

				uint32 n_threads_for_small_bins_running;
				auto n_threads_for_small_bins = n_threads - n_threads_for_big_bins;

				using SORTER_T = CRadixSorterMSD<RECORD_T, COUNTER_TYPE>;
				std::vector<std::unique_ptr<SORTER_T>> sorters;
				auto use_queue_min_recs = n_recs / 4096;
				for (n_threads_for_small_bins_running = 0; n_threads_for_small_bins_running < n_threads_for_small_bins; ++n_threads_for_small_bins_running)
				{
					sorters.emplace_back(std::make_unique<SORTER_T>(tasks_queue, 
						use_queue_min_recs, buffers[n_threads_for_small_bins_running]));
					threads.emplace_back(std::ref(*sorters.back().get()));
				}
				
				for (auto& big_bin : big_bins)
					Sort<COUNTER_TYPE>(std::get<0>(big_bin), std::get<1>(big_bin), std::get<2>(big_bin), byte - 1, n_threads_for_big_bins,
						false, is_big_threshold, n_total_recs);
				
				for (; n_threads_for_small_bins_running < n_threads; ++n_threads_for_small_bins_running)
				{
					sorters.emplace_back(std::make_unique<SORTER_T>(tasks_queue, 
						use_queue_min_recs, buffers[n_threads_for_small_bins_running]));
					threads.emplace_back(std::ref(*sorters.back().get()));
				}

				for (auto& th : threads)
					th.join();				
			}
		}
	};

	template<typename RECORD_T>
	void RadixSortMSD_template(RECORD_T* data, RECORD_T* tmp, uint64_t n_recs, uint32_t rec_size_in_bytes, uint32_t key_size_in_bytes, uint32_t n_threads)
	{
		uint32 byte = key_size_in_bytes - 1;		
		CRaduls<RECORD_T> raduls(n_threads);

		uint64 is_big_threshold = 2 * n_recs / (3 * n_threads);

		if (is_big_threshold < 4 * n_recs / 256)
			is_big_threshold = 4 * n_recs / 256;

		if (n_recs >= (1ull << 32))			
			raduls.template Sort<uint64>(data, tmp, n_recs, byte, n_threads, true, is_big_threshold, n_recs);
		else			
			raduls.template Sort<uint32>(data, tmp, n_recs, byte, n_threads, true, is_big_threshold, n_recs);
	}

	struct SortParams
	{
		uint8_t* input, *tmp;
		uint64_t n_recs;
		uint32_t key_size, rec_size, n_threads;
	};

	template<uint32_t REC_SIZE_IN_UINT64> class RecSizeDispatcher;

	template<uint32_t REC_SIZE_IN_UINT64, uint32_t KEY_SIZE_IN_UINT64>
	class KeySizeDispatcher
	{
		friend class RecSizeDispatcher<REC_SIZE_IN_UINT64>;
		friend class KeySizeDispatcher<REC_SIZE_IN_UINT64, KEY_SIZE_IN_UINT64 + 1>;
		static void Dispatch(const SortParams& p)
		{
			if ((p.key_size + 7) / 8 == KEY_SIZE_IN_UINT64)
			{
				using record_type = Record<REC_SIZE_IN_UINT64, KEY_SIZE_IN_UINT64>;
				record_type* data = reinterpret_cast<record_type*>(p.input);
				record_type* tmp = reinterpret_cast<record_type*>(p.tmp);
				RadixSortMSD_template(data, tmp, p.n_recs, p.rec_size, p.key_size, p.n_threads);
			}
			else
				KeySizeDispatcher<REC_SIZE_IN_UINT64, KEY_SIZE_IN_UINT64 - 1>::Dispatch(p);
		}
	};

	template<uint32_t REC_SIZE_IN_UINT64>
	class KeySizeDispatcher<REC_SIZE_IN_UINT64, 1>
	{
		friend class RecSizeDispatcher<REC_SIZE_IN_UINT64>;
		friend class KeySizeDispatcher<REC_SIZE_IN_UINT64, 2>;
		static void Dispatch(const SortParams& p)
		{
			using record_type = Record<REC_SIZE_IN_UINT64, 1>;
			record_type* data = reinterpret_cast<record_type*>(p.input);
			record_type* tmp = reinterpret_cast<record_type*>(p.tmp);
			RadixSortMSD_template(data, tmp, p.n_recs, p.rec_size, p.key_size, p.n_threads);
		}
	};

	template<uint32_t REC_SIZE_IN_UINT64>
	class RecSizeDispatcher
	{
		friend void RadixSortMSD(uint8_t* input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_threads);
		friend class RecSizeDispatcher<REC_SIZE_IN_UINT64 + 1>;
		static void Dispatch(const SortParams& p)
		{
			if (p.rec_size / 8 == REC_SIZE_IN_UINT64)
				KeySizeDispatcher<REC_SIZE_IN_UINT64, REC_SIZE_IN_UINT64>::Dispatch(p);
			else
				RecSizeDispatcher<REC_SIZE_IN_UINT64 - 1>::Dispatch(p);
		}
	};

	template<>
	class RecSizeDispatcher<1>
	{
		friend void RadixSortMSD(uint8_t* input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_threads);
		friend class RecSizeDispatcher<2>;
		static void Dispatch(const SortParams& p)
		{		
			KeySizeDispatcher<1, 1>::Dispatch(p);	
		}
	};

	//Non template wrapper
	//input and tmp must be aligned
	void RadixSortMSD(uint8_t* input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_threads)
	{ 
		//asserts
		if (reinterpret_cast<std::uintptr_t>(input) % ALIGNMENT)
			throw exceptions::InputNotAlignedException(ALIGNMENT);
		if (reinterpret_cast<std::uintptr_t>(tmp) % ALIGNMENT)
			throw exceptions::TempNotAlignedException(ALIGNMENT);
		if (rec_size % 8)
			throw exceptions::RecSizeNotMultipleOf8Exception();
		if (key_size > rec_size)
			throw exceptions::KeySizeGreaterThanRecSizeException();
		if (rec_size > MAX_REC_SIZE_IN_BYTES)
			throw exceptions::UsupportedRecSizeException();


		//let's go
		SortParams p;
		p.input = input;
		p.tmp = tmp;
		p.n_recs = n_recs;
		p.rec_size = rec_size;
		p.key_size = key_size;
		p.n_threads = n_threads;

		RecSizeDispatcher<MAX_REC_SIZE_IN_BYTES / 8>::Dispatch(p);		
	}

	void CleanTmpArray(uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t n_threads)
	{
		std::vector<std::thread> cleaning_threads;
		auto n_bytes = rec_size * n_recs;
		auto part_size = n_bytes / n_threads;
		for (uint32_t th_id = 0; th_id < n_threads; ++th_id)
		{
			cleaning_threads.emplace_back([th_id, n_bytes, part_size, tmp]
			{
				auto start = th_id * part_size;
				auto end = start + part_size;
				if (end > n_bytes)
					end = n_bytes;
				for (uint64_t i = start; i < end; i += 4096)
					tmp[i] = 0;
			});
		}
		for (auto& t : cleaning_threads)
			t.join();
	}
}

