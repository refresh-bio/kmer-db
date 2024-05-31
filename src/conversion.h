#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek
*/

#include <memory>
#include <cassert>
#include <string>
#include <cstdint>
#include <cstring>
#include <type_traits>

// ************************************************************************************
class NumericConversions
{
public:
	static char digits[100000 * 5];
	static uint64_t powers10[15];
	struct _si {
		_si()
		{
			for (int i = 0; i < 100000; ++i)
			{
				int dig = i;

				digits[i * 5 + 4] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 3] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 2] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 1] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 0] = '0' + dig;
			}

			powers10[0] = 1;
			for (int i = 1; i < 15; ++i)
				powers10[i] = 10 * powers10[i - 1];
		}
	} static _init;

	static int NDigits(uint64_t v)
	{
		return (v < 10000)
			? (v < 100 ? (v < 10 ? 1 : 2) : (v < 1000 ? 3 : 4))
			: (v < 1000000 ? (v < 100000 ? 5 : 6) : (v < 10000000 ? 7 : 8));
	}

	// Works only for len <= 5
	static void short_str_cpy_upto_5(char* dest, char* src, size_t len)
	{
		assert(len <= 5);

		if (len == 1)
			dest[0] = src[0];
		else if (len == 2)
		{
			dest[0] = src[0];
			dest[1] = src[1];
		}
		else if (len == 3)
		{
			dest[0] = src[0];
			dest[1] = src[1];
			dest[2] = src[2];
		}
		else if (len == 4)
		{
			dest[0] = src[0];
			dest[1] = src[1];
			dest[2] = src[2];
			dest[3] = src[3];
		}
		else
		{
			dest[0] = src[0];
			dest[1] = src[1];
			dest[2] = src[2];
			dest[3] = src[3];
			dest[4] = src[4];
		}

		return;
	}

	static void short_str_cpy_5(char* dest, char* src)
	{
		dest[0] = src[0];
		dest[1] = src[1];
		dest[2] = src[2];
		dest[3] = src[3];
		dest[4] = src[4];
	}

	static int Int2PChar(uint64_t val, char *str)
	{
		if (val >= 1000000000000000ull)
		{
			uint64_t dig1 = val / 1000000000000000ull;
			val -= dig1 * 1000000000000000ull;
			uint64_t dig2 = val / 10000000000ull;
			val -= dig2 * 10000000000ull;
			uint64_t dig3 = val / 100000ull;
			uint64_t dig4 = val - dig3 * 100000ull;

			int ndig = NDigits(dig1);

/*			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);
			std::memcpy(str + ndig + 10, digits + dig4 * 5, 5);*/
			short_str_cpy_upto_5(str, digits + dig1 * 5 + (5 - ndig), ndig);
			short_str_cpy_5(str + ndig, digits + dig2 * 5);
			short_str_cpy_5(str + ndig + 5, digits + dig3 * 5);
			short_str_cpy_5(str + ndig + 10, digits + dig4 * 5);

			return ndig + 15;
		}
		else if (val >= 10000000000ull)
		{
			uint64_t dig1 = val / 10000000000ull;
			val -= dig1 * 10000000000ull;
			uint64_t dig2 = val / 100000ull;
			uint64_t dig3 = val - dig2 * 100000ull;

			int ndig = NDigits(dig1);

/*			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);*/
			short_str_cpy_upto_5(str, digits + dig1 * 5 + (5 - ndig), ndig);
			short_str_cpy_5(str + ndig, digits + dig2 * 5);
			short_str_cpy_5(str + ndig + 5, digits + dig3 * 5);

			return ndig + 10;
		}
		else if (val >= 100000ull)
		{
			uint64_t dig1 = val / 100000ull;
			uint64_t dig2 = val - dig1 * 100000ull;

			int ndig = NDigits(dig1);

//			memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
//			memcpy(str + ndig, digits + dig2 * 5, 5);
			short_str_cpy_upto_5(str, digits + dig1 * 5 + (5 - ndig), ndig);
//			short_str_cpy(str + ndig, digits + dig2 * 5, 5);
			short_str_cpy_5(str + ndig, digits + dig2 * 5);

			return ndig + 5;
		}
		else
		{
			int ndig = NDigits(val);

//			memcpy(str, digits + val * 5 + (5 - ndig), ndig);
			short_str_cpy_upto_5(str, digits + val * 5 + (5 - ndig), ndig);

			return ndig;
		}
	}

	static int Double2PChar(double val, uint32_t prec, char *str)
	{
/*		int64_t a = (int64_t)val;
		int64_t b = (int64_t)((1.0 + (val - (double)a)) * powers10[prec] + 0.5);

		int r1 = Int2PChar(a, str);
		int r2 = Int2PChar(b, str + r1);
		str[r1] = '.';

		return r1 + r2;*/

		int neg = 0;

		if (val < 0)
		{
			*str++ = '-';
			val = -val;
			neg = 1;
		}

		uint64_t x = (uint64_t) (val * powers10[prec] + 0.5);

		if (x < powers10[prec])						// |val| < 1.0
		{
			*str++ = '0';
			*str++ = '.';

			int to_move = Int2PChar(x, str);
			int shift = prec - to_move;

			if (shift)
			{
				for (int i = 0; i < to_move; ++i)
					str[prec - i - 1] = str[prec - i - 1 - shift];
				for (int i = 0; i < shift; ++i)
					str[i] = '0';
			}

			return prec + 2 + neg;
		}
		else
		{
			int r = Int2PChar(x, str + 1);

			int to_move = r - prec;

			for (int i = 0; i < to_move; ++i)
				str[i] = str[i+1];
			str[to_move] = '.';

			return r + 1 + neg;
		}
	}

	static long int strtol(const char* str, char** endptr) {
		long int val = 0;
		char* p = (char*)str;
		bool is_negative = false;


		if (*p == '-')
		{
			is_negative = true;
			++p;
		}

		while (*p >= '0' && *p <= '9')
		{
			val = val * 10 + (*p++ - '0');
		}

		if (endptr)
			*endptr = p;

		return is_negative ? -val : val;
	}
};


// integral specialization
template <typename Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type* = nullptr>
int num2str(Integer val, char *out) {
	return NumericConversions::Int2PChar((uint64_t)val, out);
}

// floating point specialization
template <typename Floating, typename std::enable_if<std::is_floating_point<Floating>::value, int>::type* = nullptr>
int num2str(Floating val, char *out) {
	if (val == 0) {
		*out = '0';
		return 1;
	}
	return NumericConversions::Double2PChar((double)val, 6, out);
}

// pair specialization
template <typename T, typename U>
int num2str(const std::pair<T,U> val, char *out) {
	char* ptr = out;
	ptr += num2str(val.first, ptr);
	*ptr++ = ':';
	ptr += num2str(val.second, ptr);
	
	return ptr - out;
}

// collection specialization
template <typename T>
int num2str(const T* collection, size_t size, char delim, char* out) {
	char* ptr = out;
	for (size_t i = 0; i < size; ++i) {
		ptr += num2str(*collection++, ptr);
		*ptr++ = delim;
	}

	return ptr - out;
}

template <typename T>
int num2str_sparse(const T* collection, size_t size, char delim, char* out, const T sparse_val = 0) {
	char* ptr = out;
	for (size_t i = 0; i < size; ++i, collection++) {
		if (*collection != sparse_val) {
			ptr += num2str(i + 1, ptr);
			*ptr++ = ':';
			ptr += num2str(*collection, ptr);
			*ptr++ = delim;
		}
	}

	return ptr - out;
}