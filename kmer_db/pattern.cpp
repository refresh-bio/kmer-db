#include "pattern.h"


CEliasGamma pattern_t::elias;


char* pattern_t::pack(char* buffer) const {
	// store members
	*reinterpret_cast<decltype(num_kmers)*>(buffer) = num_kmers;
	buffer += sizeof(decltype(num_kmers));

	*reinterpret_cast<decltype(parent_id)*>(buffer) = parent_id;
	buffer += sizeof(decltype(parent_id));

	*reinterpret_cast<decltype(is_parent)*>(buffer) = is_parent;
	buffer += sizeof(decltype(is_parent));

	*reinterpret_cast<decltype(num_samples)*>(buffer) = num_samples;
	buffer += sizeof(decltype(num_samples));

	*reinterpret_cast<decltype(num_local_samples)*>(buffer) = num_local_samples;
	buffer += sizeof(decltype(num_local_samples));

	*reinterpret_cast<decltype(num_bits)*>(buffer) = num_bits;
	buffer += sizeof(decltype(num_bits));

	*reinterpret_cast<decltype(last_sample_id)*>(buffer) = last_sample_id;
	buffer += sizeof(decltype(last_sample_id));

	size_t data_bytes = get_data_bytes();
	if (data_bytes) {
		memcpy(buffer, reinterpret_cast<char*>(data), data_bytes);
		buffer += data_bytes;
	}

	return buffer;
}

char * pattern_t::unpack(char* buffer) {
	if (num_local_samples) {
		delete[] data;
	}

	num_kmers = *reinterpret_cast<decltype(num_kmers)*>(buffer);
	buffer += sizeof(decltype(num_kmers));

	parent_id = *reinterpret_cast<decltype(parent_id)*>(buffer);
	buffer += sizeof(decltype(parent_id));

	is_parent = *reinterpret_cast<decltype(is_parent)*>(buffer);
	buffer += sizeof(decltype(is_parent));

	num_samples = *reinterpret_cast<decltype(num_samples)*>(buffer);
	buffer += sizeof(decltype(num_samples));

	num_local_samples = *reinterpret_cast<decltype(num_local_samples)*>(buffer);
	buffer += sizeof(decltype(num_local_samples));

	num_bits = *reinterpret_cast<decltype(num_bits)*>(buffer);
	buffer += sizeof(decltype(num_bits));

	last_sample_id = *reinterpret_cast<decltype(last_sample_id)*>(buffer);
	buffer += sizeof(decltype(last_sample_id));


	size_t num_bytes = get_data_bytes();

	if (num_bytes) {
		data = new uint64_t[num_bytes / sizeof(uint64_t)];
		std::memcpy(reinterpret_cast<char*>(data), buffer, num_bytes);
		buffer += num_bytes;
	}

	return buffer;
}

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