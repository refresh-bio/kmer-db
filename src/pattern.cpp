/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "pattern.h"

CEliasGamma pattern_t::elias;

// *****************************************************************************************
//
char* pattern_t::pack(char* buffer) const {
	// store members
	memcpy(buffer, &num_kmers, sizeof(int64_t));
	buffer += sizeof(int64_t);

	memcpy(buffer, &parent_id, sizeof(int64_t));
	buffer += sizeof(int64_t);

	memcpy(buffer, &num_samples, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);

	memcpy(buffer, &num_local_samples, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);

	memcpy(buffer, &last_sample_id, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);
	
	memcpy(buffer, &num_bits, sizeof(uint32_t));
	buffer += sizeof(uint32_t);

	uint64_t tmp = (uint64_t)is_parent;
	memcpy(buffer, &tmp, sizeof(uint32_t));
	buffer += sizeof(uint64_t);

	size_t data_bytes = get_data_bytes();
	if (data_bytes) {
		memcpy(buffer, reinterpret_cast<char*>(data), data_bytes);
		buffer += data_bytes;
	}

	return buffer;
}

// *****************************************************************************************
//
char * pattern_t::unpack(char* buffer) {
	if (num_local_samples) {
#ifdef USE_MALLOC
		free(data);
#else
		delete[] data;
#endif
	}

	memcpy(&num_kmers, buffer, sizeof(int64_t));
	buffer += sizeof(int64_t);

	memcpy(&parent_id, buffer, sizeof(int64_t));
	buffer += sizeof(int64_t);

	memcpy(&num_samples, buffer, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);

	memcpy(&num_local_samples, buffer, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);

	memcpy(&last_sample_id, buffer, sizeof(sample_id_t));
	buffer += sizeof(sample_id_t);

	memcpy(&num_bits, buffer, sizeof(uint32_t));
	buffer += sizeof(uint32_t);

	uint64_t tmp;
	memcpy(&tmp, buffer, sizeof(uint32_t));
	is_parent = tmp;
	buffer += sizeof(uint64_t);

	size_t num_bytes = get_data_bytes();

	if (num_bytes) {
#ifdef USE_MALLOC
		data = (uint64_t*)malloc(num_bytes);
#else
		data = new uint64_t[num_bytes / sizeof(uint64_t)];
#endif
		std::memcpy(reinterpret_cast<char*>(data), buffer, num_bytes);
		buffer += num_bytes;
	}

	return buffer;
}

// *****************************************************************************************
//
void pattern_t::decodeSamples(uint32_t* out) const {
	if (num_local_samples) {
		out[num_local_samples - 1] = last_sample_id;
		if (num_local_samples > 1) {
			elias.Decode(data, num_bits, out);
			for (int i = num_local_samples - 2; i >= 0; --i) {
				out[i] = out[i + 1] - out[i];
			}
		}
	}
}