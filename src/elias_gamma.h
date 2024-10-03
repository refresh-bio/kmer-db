#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <iostream>

using namespace std;

#ifdef __GNUG__
#define FORCE_INLINE inline __attribute__((always_inline))
using __int64 = long long int;
#else
#define FORCE_INLINE __forceinline
#endif 

#define USE_MALLOC

// *****************************************************************************************
//
class CEliasGamma
{
	const uint32_t LEN_ILOG2 = 15;
	const uint32_t LEN_PREFIX_LEN = 15;

	uint64_t shifts[64];
	uint64_t masks[65];
	uint8_t *lut_ilog2;
	uint8_t *lut_prefix_len;

	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t round16(uint32_t x)
	{
		return (x + 15) / 16 * 16;
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t ilog2(uint32_t x)
	{
		uint32_t r;

		if (x < (1u << LEN_ILOG2))
			return lut_ilog2[x];
		else
		{
			x >>= LEN_ILOG2;
			r = LEN_ILOG2;
		}

		for(; x >= 64; r += 6)
			x >>= 6;
		
		for(; x; ++r)
			x >>= 1;
		
		return r;
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t no_of_16B_blocks(uint32_t no_bits)
	{
		if(no_bits <= 256)
			return (no_bits + 127) / 128;
		else if(no_bits <= 512)
			return 2 * ((no_bits + 255) / 256);
		else
			return 3 * ((no_bits + 127 + 256) / (128 + 256));
	}

	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t calculate_size(uint32_t value)
	{
		return 2 * ilog2(value) - 1;
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t calculate_size(uint32_t *input, uint32_t input_size)
	{
		uint32_t r = 0;
		
		for(uint32_t i = 0; i < input_size; ++i)
			r += calculate_size(input[i]);
		
		return r;
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	void encode(uint32_t value, uint64_t *&output, uint32_t &output_size_in_bits)
	{
		uint32_t prefix_len = ilog2(value);
		uint32_t code_len = 2 * prefix_len - 1;
		
		uint64_t code = masks[prefix_len - 1];				// necessary no. of 1s
		code <<= prefix_len;								// 0 bit and space for remaining bits
		code += value - shifts[prefix_len - 1];				// remaining bits

		uint32_t current_block_filled_bits = output_size_in_bits % 64;
		uint32_t current_block_word_id = output_size_in_bits / 64;
		
		uint32_t no_empty_bits = 64 - current_block_filled_bits;
		
		if(code_len <= no_empty_bits)
			output[current_block_word_id] += code << (no_empty_bits - code_len);
		else
		{
			uint32_t no_remaining_bits = code_len - no_empty_bits;
			output[current_block_word_id] += code >> no_remaining_bits;
			output[current_block_word_id + 1] = (code & masks[no_remaining_bits]) << (64 - no_remaining_bits);
		}

		output_size_in_bits += code_len;
	}

	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t decode(uint64_t *input, uint32_t &input_pos_in_bits)
	{
#ifdef NAIVE_DECODE
		while (true)
		{
			uint32_t word_id = input_pos_in_bits / 64;
			uint32_t bit_id = 63 - input_pos_in_bits % 64;
			++input_pos_in_bits;

			if (input[word_id] & shifts[bit_id])
				++prefix_len;
			else
				break;
		}

		uint32_t r = 1;
		for (uint32_t i = 0; i < prefix_len; ++i)
		{
			uint32_t word_id = input_pos_in_bits / 64;
			uint32_t bit_id = 63 - input_pos_in_bits % 64;

			r = (r << 1) + ((input[word_id] >> bit_id) & 1);

			++input_pos_in_bits;
		}
#else
		uint32_t word_id = input_pos_in_bits / 64;
		uint32_t bit_id = 63 - input_pos_in_bits % 64;

		// 1 is frequent, so quick test for it
		if (!(input[word_id] & shifts[bit_id]))
		{
			++input_pos_in_bits;
			return 1;
		}

		uint32_t prefix_len = 0;
		bool prefix_ready = false;

		while(!prefix_ready)
			if (bit_id + 1 >= LEN_PREFIX_LEN)		
			{
				uint64_t lut_id = input[word_id] >> (bit_id - LEN_PREFIX_LEN + 1);
				lut_id &= masks[LEN_PREFIX_LEN];

				uint32_t iter_prefix_len = lut_prefix_len[lut_id];
				if (iter_prefix_len < LEN_PREFIX_LEN)
				{
					prefix_len += iter_prefix_len;
					input_pos_in_bits += prefix_len + 1;
					bit_id -= iter_prefix_len + 1;
					prefix_ready = true;
				}
				else
				{
					bit_id -= iter_prefix_len;
					prefix_len += iter_prefix_len;
				}
			}
			else 
			{
				uint64_t lut_id = input[word_id] & masks[bit_id + 1];
				lut_id <<= LEN_PREFIX_LEN - bit_id - 1;

				uint32_t iter_prefix_len = lut_prefix_len[lut_id];
				if (iter_prefix_len < bit_id + 1)
				{
					prefix_len += iter_prefix_len;
					input_pos_in_bits += prefix_len + 1;
					bit_id -= iter_prefix_len + 1;
					prefix_ready = true;
				}
				else
				{
					bit_id -= iter_prefix_len;
					prefix_len += iter_prefix_len;
				}
				break;
			}

		if (!prefix_ready)	
		{
			++word_id;
			bit_id = 63;

			while (!prefix_ready)
			{
				uint64_t lut_id = input[word_id] >> (bit_id - LEN_PREFIX_LEN + 1);
				lut_id &= masks[LEN_PREFIX_LEN];

				uint32_t iter_prefix_len = lut_prefix_len[lut_id];
				if (iter_prefix_len < LEN_PREFIX_LEN)
				{
					prefix_len += iter_prefix_len;
					input_pos_in_bits += prefix_len + 1;
					bit_id -= iter_prefix_len + 1;
					prefix_ready = true;
				}
				else
				{
					prefix_len += iter_prefix_len;
					bit_id -= iter_prefix_len;
				}
			}
		}

		uint32_t r;

		if (bit_id + 1 >= prefix_len)
		{
			input_pos_in_bits += prefix_len;
			r = (uint32_t) (shifts[prefix_len] + ((input[word_id] >> (bit_id - prefix_len + 1)) & masks[prefix_len]));
		}
		else
		{
			input_pos_in_bits += prefix_len;
			r = (uint32_t) (shifts[prefix_len]);
			prefix_len -= bit_id + 1;
			r += (uint32_t) ((input[word_id] & masks[bit_id + 1]) << prefix_len);
			r += (uint32_t) (input[word_id + 1] >> (64 - prefix_len));
		}
#endif
		return r;
	}

public:
	// *****************************************************************************************
	//
	CEliasGamma()
	{
		//  1 << i
		shifts[0] = 1;
		for (int i = 1; i < 64; ++i)
			shifts[i] = shifts[i - 1] << 1;

		// masks[i] contains i set bits
		masks[0] = 0;
		for (int i = 1; i < 65; ++i)
			masks[i] = (masks[i - 1] << 1) + 1ull;

		// lut_ilog2[i] = min. no. of bits for number i
		lut_ilog2 = new uint8_t[1ull << LEN_ILOG2];
		lut_ilog2[0] = 0;
		lut_ilog2[1] = 1;

		for (int i = 2; i < (1 << LEN_ILOG2); ++i)
			lut_ilog2[i] = lut_ilog2[i / 2] + 1;

		// lut_prefix_len[i] - no. of leading 1 in LEN_PREFIX_LEN lower bits of i
		lut_prefix_len = new uint8_t[1ull << LEN_PREFIX_LEN];

		uint32_t prefix_len = 0;
		uint32_t pos = 0;
		for (uint32_t range = 1 << (LEN_PREFIX_LEN - 1); range; range >>= 1)
		{
			for (uint32_t i = 0; i < range; ++i)
				lut_prefix_len[pos + i] = prefix_len;

			pos += range;
			++prefix_len;
		}

		lut_prefix_len[pos] = LEN_PREFIX_LEN;
	}

	// *****************************************************************************************
	//
	~CEliasGamma()
	{
		delete[] lut_ilog2;
		delete[] lut_prefix_len;
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	void Encode(uint32_t *input, uint32_t input_size, uint64_t *&output, uint32_t &output_size_in_bits)
	{
		uint32_t predicted_size_in_bits = calculate_size(input, input_size);
		uint32_t predicted_size_in_8B = no_of_16B_blocks(predicted_size_in_bits) * 2;

#ifdef USE_MALLOC
		if (output)
			free(output);
		output = (uint64_t*)calloc(predicted_size_in_8B, sizeof(uint64_t));
#else
		if (output)
			delete[] output;
		output = new uint64_t[predicted_size_in_8B];
		fill_n(output, predicted_size_in_8B, 0);
#endif

		output_size_in_bits = 0;
		for (uint32_t i = 0; i < input_size; ++i)
			encode(input[i], output, output_size_in_bits);
	}

	// *****************************************************************************************
	//
	FORCE_INLINE
	void Encode(uint32_t value, uint64_t *&output, uint32_t &output_size_in_bits)
	{
		uint32_t size = calculate_size(value);
		uint32_t predicted_size_in_bits = size + output_size_in_bits;
		uint32_t current_no_blocks = no_of_16B_blocks(output_size_in_bits);
		uint32_t new_no_blocks = no_of_16B_blocks(predicted_size_in_bits);
		
		// Do we need to reallocate?
		if(current_no_blocks != new_no_blocks)
		{
#ifdef USE_MALLOC
			output = (uint64_t*)realloc(output, new_no_blocks * 16);
			fill_n(output + current_no_blocks * 2, new_no_blocks * 2 - current_no_blocks * 2, 0);
#else
			uint64_t *old_output = output;
			output = new uint64_t[new_no_blocks * 2];
			memcpy(output, old_output, current_no_blocks * 16);
			fill_n(output + current_no_blocks * 2, new_no_blocks * 2 - current_no_blocks * 2, 0);

			delete[] old_output;
#endif
		}
		
		encode(value, output, output_size_in_bits);
	}
	
	// *****************************************************************************************
	//
	FORCE_INLINE
	uint32_t Decode(uint64_t *input, uint32_t &input_pos_in_bits)
	{
		return decode(input, input_pos_in_bits);
	}
	
	// *****************************************************************************************
	//
	// Assumption: memory is allocated
	FORCE_INLINE
	void Decode(uint64_t *input, uint32_t input_size_in_bits, uint32_t *output)
	{
		uint32_t input_pos_in_bits = 0;
		uint32_t i = 0;
		
		while (input_pos_in_bits < input_size_in_bits)
			output[i++] = decode(input, input_pos_in_bits);
	}
};
