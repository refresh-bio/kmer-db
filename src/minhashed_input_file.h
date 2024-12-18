#pragma once
#include "input_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "kmer_extract.h"
#include "parallel_sorter.h"

#include "../libs/refresh/compression/lib/file_wrapper.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>
#include <filesystem>

#ifdef USE_RADULS
#include <raduls.h>
#endif

// *****************************************************************************************
//
class MihashedInputFile : public InputFile {
public:
	inline bool open(const std::string& filename) override;

	inline bool load(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue) override;

	inline bool store(
		const std::string& filename, 
		const kmer_t* kmers, 
		size_t kmersCount, 
		uint32_t kmerLength, 
		double filterValue);

protected:
	const uint32_t MINHASH_FORMAT_SIGNATURE = 0xfedcba98;

	std::vector<kmer_t> kmers;

	uint32_t kmerLength{ 0 };

	double fraction{ 0 };

	bool status{ false };
};



// *****************************************************************************************
//
bool MihashedInputFile::open(const std::string& filename) {
	std::ifstream file(filename + ".minhash", std::ios_base::binary);
	status = false;
	if (file) {
		uint32_t signature = 0;
		file.read(reinterpret_cast<char*>(&signature), sizeof(uint32_t));
		if (signature == MINHASH_FORMAT_SIGNATURE) {
			size_t numKmers;
			file.read(reinterpret_cast<char*>(&numKmers), sizeof(size_t));
			kmers.resize(numKmers);

			file.read(reinterpret_cast<char*>(kmers.data()), sizeof(kmer_t) * numKmers);
			file.read(reinterpret_cast<char*>(&kmerLength), sizeof(kmerLength));
			file.read(reinterpret_cast<char*>(&fraction), sizeof(fraction));
			if (file) {
				status = true;
			}
		}
		file.close();
	}

	if (!status) {
		kmers.clear();
	}

	return status;
}

// *****************************************************************************************
//
bool MihashedInputFile::load(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue) {
	if (!status) {
		return false;
	}

	kmersBuffer = std::move(this->kmers);
	kmers = kmersBuffer.data();
	kmersCount = kmersBuffer.size();
	kmerLength = this->kmerLength;
	filterValue = this->fraction;
	return true;
}

// *****************************************************************************************
//
bool MihashedInputFile::store(const std::string& filename, const kmer_t* kmers, size_t kmersCount, uint32_t kmerLength, double filterValue) {
	ofstream ofs(filename + ".minhash", std::ios_base::binary);
	ofs.write(reinterpret_cast<const char*>(&MINHASH_FORMAT_SIGNATURE), sizeof(MINHASH_FORMAT_SIGNATURE));

	ofs.write(reinterpret_cast<const char*>(&kmersCount), sizeof(size_t));
	ofs.write(reinterpret_cast<const char*>(kmers), kmersCount * sizeof(kmer_t));
	ofs.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength));
	ofs.write(reinterpret_cast<const char*>(&filterValue), sizeof(filterValue));
	return true;
}



