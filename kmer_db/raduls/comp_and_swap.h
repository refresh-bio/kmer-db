#pragma once
#include "defs.h"
#include "xmmintrin.h"
namespace raduls
{
	namespace small_sort
	{
		template<typename RECORD_T, typename Comp>
		struct SwapLowerGreater
		{
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				const auto a = Comp{}(lhs, rhs) ? lhs : rhs;
				const auto b = Comp{}(rhs, lhs) ? lhs : rhs;
				lhs = a;
				rhs = b;
			}
		};

		template<unsigned RECORD_SIZE, typename Comp>
		struct SwapLowerGreater<Record<RECORD_SIZE, 2>, Comp>
		{
			using RECORD_T = Record<RECORD_SIZE, 2>;
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				const auto a = Comp{}(lhs, rhs) ? lhs : rhs;
				const auto b = Comp{}(rhs, lhs) ? lhs : rhs;
				lhs = a;
				rhs = b;
			}
		};

		template<unsigned RECORD_SIZE, typename Comp>
		struct SwapLowerGreater<Record<RECORD_SIZE, 1>, Comp>
		{
			using RECORD_T = Record<RECORD_SIZE, 1>;
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				const auto a = lhs.data[0] < rhs.data[0] ? lhs : rhs;
				const auto b = lhs.data[0] > rhs.data[0] ? lhs : rhs;
				lhs = a;
				rhs = b;
			}
		};

		template<typename Comp>
		struct SwapLowerGreater<Record<2, 1>, Comp>
		{
			using RECORD_T = Record<2, 1>;
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				__m128i a = _mm_loadu_si128(Comp{}(lhs, rhs) ? (const __m128i*)lhs.data : (const __m128i*)rhs.data);
				__m128i b = _mm_loadu_si128(!Comp{}(lhs, rhs) ? (const __m128i*)lhs.data : (const __m128i*)rhs.data);
				_mm_storeu_si128((__m128i*)lhs.data, a);
				_mm_storeu_si128((__m128i*)rhs.data, b);
			}
		};

		template<typename Comp>
		struct SwapLowerGreater<Record<1, 1>, Comp>
		{
			using RECORD_T = Record<1, 1>;
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				const auto a = lhs.data < rhs.data ? lhs.data : rhs.data;
				const auto b = lhs.data > rhs.data ? lhs.data : rhs.data;
				lhs.data = a;
				rhs.data = b;
			}
		};

		template<typename RECORD_T>
		struct LessFirstNotEqual
		{
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				for (int i = RECORD_T::KEY_SIZE - 1; i >= 0; --i)
					if (lhs.data[i] != rhs.data[i])
						return lhs.data[i] < rhs.data[i];
				return false;
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstNotEqual<Record<RECORD_SIZE, 3>>
		{
			using RECORD_T = Record<RECORD_SIZE, 3>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				if (lhs.data[2] != rhs.data[2])
					return lhs.data[2] < rhs.data[2];
				if (lhs.data[1] != rhs.data[1])
					return lhs.data[1] < rhs.data[1];
				return lhs.data[0] < rhs.data[0];
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstNotEqual<Record<RECORD_SIZE, 2>>
		{
			using RECORD_T = Record<RECORD_SIZE, 2>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data[1] == rhs.data[1] ? lhs.data[0] < rhs.data[0] : lhs.data[1] < rhs.data[1];
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstNotEqual<Record<RECORD_SIZE, 1>>
		{
			using RECORD_T = Record<RECORD_SIZE, 1>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data[0] < rhs.data[0];
			}
		};

		template<>
		struct LessFirstNotEqual<Record<1, 1>>
		{
			using RECORD_T = Record<1, 1>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data < rhs.data;
			}
		};

		template<typename RECORD_T>
		struct LessFirstLower
		{
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				for (int32 i = RECORD_T::KEY_SIZE - 1; i >= 0; --i)
					if (lhs.data[i] < rhs.data[i])
						return true;
					else if (lhs.data[i] > rhs.data[i])
						return false;
				return false;
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstLower<Record<RECORD_SIZE, 3>>
		{
			using RECORD_T = Record<RECORD_SIZE, 3>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				if (lhs.data[2] < rhs.data[2])
					return true;
				if (lhs.data[2] == rhs.data[2])
				{
					if (lhs.data[1] < rhs.data[1])
						return true;
					if (lhs.data[1] == rhs.data[1])
						return lhs.data[0] < rhs.data[0];
				}
				return false;
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstLower<Record<RECORD_SIZE, 2>>
		{
			using RECORD_T = Record<RECORD_SIZE, 2>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data[1] < rhs.data[1] ? true : lhs.data[1] == rhs.data[1] ? lhs.data[0] < rhs.data[0] : false;
			}
		};

		template<unsigned RECORD_SIZE>
		struct LessFirstLower<Record<RECORD_SIZE, 1>>
		{
			using RECORD_T = Record<RECORD_SIZE, 1>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data[0] < rhs.data[0];
			}
		};

		template<>
		struct LessFirstLower<Record<1, 1>>
		{
			using RECORD_T = Record<1, 1>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data < rhs.data;
			}
		};

		template<typename RECORD_T>
		struct FinishingWrapper
		{
			FORCE_INLINE void operator()(RECORD_T* data, uint32 size)
			{
				ins_sort1a(data, data + size, LessFirstLower<RECORD_T>{});
			}
		};

		template<typename RECORD_T>
		struct FakeFinishingSorter
		{
			FORCE_INLINE void operator()(RECORD_T*, uint32) {};
		};

		template<typename RECORD_T>
		struct MS_uint64_lower
		{
			bool operator()(const RECORD_T& lhs, const RECORD_T& rhs);
		};

		template<unsigned RECORD_SIZE, unsigned KEY_SIZE>
		struct MS_uint64_lower<Record<RECORD_SIZE, KEY_SIZE>>
		{
			using RECORD_T = Record<RECORD_SIZE, KEY_SIZE>;
			FORCE_INLINE bool operator()(const RECORD_T& lhs, const RECORD_T& rhs)
			{
				return lhs.data[KEY_SIZE - 1] < rhs.data[KEY_SIZE - 1];
			}
		};

		template<typename RECORD_T, typename Comp>
		struct IntrSwapper
		{
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs);
		};

		template <typename Comp>
		struct IntrSwapper<Record<2, 2>, Comp>
		{
			using RECORD_T = Record<2, 2>;
			FORCE_INLINE void operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				__m128i a = _mm_loadu_si128(Comp{}(lhs, rhs) ? (const __m128i*)lhs.data : (const __m128i*)rhs.data);
				__m128i b = _mm_loadu_si128(!Comp{}(lhs, rhs) ? (const __m128i*)lhs.data : (const __m128i*)rhs.data);
				_mm_storeu_si128((__m128i*)lhs.data, a);
				_mm_storeu_si128((__m128i*)rhs.data, b);
			}
		};

		template<typename RECORD_T>
		struct LS_uint64_lower
		{
			bool operator()(RECORD_T& lhs, RECORD_T& rhs);
		};

		template<>
		struct LS_uint64_lower<Record<2, 2>>
		{
			using RECORD_T = Record<2, 2>;
			FORCE_INLINE bool operator()(RECORD_T& lhs, RECORD_T& rhs)
			{
				return lhs.data[0] < rhs.data[0];
			}
		};
	}
}
