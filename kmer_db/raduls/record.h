#pragma once
#include "defs.h"

namespace raduls
{
	template<unsigned _RECORD_SIZE, unsigned _KEY_SIZE>
	struct Record
	{
		static const uint32 RECORD_SIZE = _RECORD_SIZE;
		static const uint32 KEY_SIZE = _KEY_SIZE;

		uint64 data[_RECORD_SIZE];

		FORCE_INLINE bool operator<(const Record<_RECORD_SIZE, _KEY_SIZE>& rhs)
		{
			for (int32 i = _KEY_SIZE - 1; i >= 0; --i)
				if (data[i] < rhs.data[i])
					return true;
				else if (data[i] > rhs.data[i])
					return false;
			return false;
		}
		
		FORCE_INLINE bool operator==(const Record<_RECORD_SIZE, _KEY_SIZE>& rhs)
		{
			for (uint32 i = 0; i < _KEY_SIZE; ++i)
				if (data[i] != rhs.data[i])
					return false;
			return true;
		}
	};

	template<unsigned _RECORD_SIZE>
	struct Record<_RECORD_SIZE, 2>
	{
		static const uint32 RECORD_SIZE = _RECORD_SIZE;
		static const uint32 KEY_SIZE = 2;

		uint64 data[_RECORD_SIZE];

		FORCE_INLINE bool operator<(const Record<_RECORD_SIZE, 2>& rhs)
		{
			return data[1] < rhs.data[1] || data[1] == rhs.data[1] && data[0] < rhs.data[0];
		}

		FORCE_INLINE bool operator==(const Record<_RECORD_SIZE, 2>& rhs) const
		{
			return data[0] == rhs.data[0] && data[1] == rhs.data[1];
		}
	};

	template<unsigned _RECORD_SIZE>
	struct Record<_RECORD_SIZE, 1>
	{
		static const uint32 RECORD_SIZE = _RECORD_SIZE;
		static const uint32 KEY_SIZE = 1;

		uint64 data[_RECORD_SIZE];

		FORCE_INLINE bool operator<(const Record<_RECORD_SIZE, 1>& rhs)
		{
			return data[0] < rhs.data[0];
		}

		FORCE_INLINE bool operator==(const Record<_RECORD_SIZE, 1>& rhs) const
		{
			return data[0] == rhs.data[0];
		}
	};

	template<>
	struct Record<1, 1>
	{
		static const uint32 RECORD_SIZE = 1;
		static const uint32 KEY_SIZE = 1;

		uint64 data;

		FORCE_INLINE bool operator<(const Record<1, 1>& rhs)
		{
			return data < rhs.data;
		}

		FORCE_INLINE bool operator==(const Record<1, 1>& x) const
		{
			return data == x.data;
		}
	};	
}