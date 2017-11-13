#pragma once
#include <cstdint>

#define RADULS_VER  "2.0.0"
#define RADULS_DATE "2017-10-22"

#ifdef __GNUG__
#define FORCE_INLINE inline __attribute__((always_inline))
using __int64 = long long int;
#else
#define FORCE_INLINE __forceinline
#endif

namespace raduls
{
	using uchar = unsigned char;
	using int32 = int32_t;
	using uint32 = uint32_t;
	using int64 = int64_t;
	using uint64 = uint64_t;
}