/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include "input_file.h"
#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "kmer_extract.h"
#include "parallel_sorter.h"

#include "../libs/refresh/file_wrapper.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>
#include <filesystem>

#ifdef USE_RADULS
	#include <raduls.h>
#endif

using namespace std;

// *****************************************************************************************
//
 bool GenomeInputFile::open(const std::string& filename) {

	vector<string> extensions { 
		"", ".fa", ".fna", ".fasta", ".fastq",
		".gz", ".fa.gz", ".fna.gz", ".fasta.gz", ".fastq.gz"};

	string extname;
	for (const auto& ext : extensions) {
		if (filesystem::exists(filename + ext)) {
			extname = filename + ext;
			break;
		}
	}

	status = false;

	refresh::stream_in_file input(extname);

	if (!input.is_open()) {
		return false;
	}

	size_t fsize = filesystem::file_size(extname);
	refresh::stream_decompression stream(&input);

	auto format = stream.get_format();
	// assume compression ratio of max 5
	size_t multiplier = (format == refresh::stream_decompression::format_t::text) ? 1 : 5;
	// !!! TODO: why allocate in parts? 
	size_t block_size = std::min(fsize * multiplier + 1, (size_t)(100ULL << 20));  // for null termination

	size_t buf_size = block_size;
	rawData = reinterpret_cast<char*>(malloc(buf_size)); 
	rawSize = 0;

	size_t n_read = 0;
	for (;;) {
		stream.read(rawData + rawSize, block_size, n_read);

		rawSize += n_read;

		// more data is coming...
		if (n_read > 0 && n_read == block_size) {
			block_size = static_cast<size_t>(block_size * 1.5);

			buf_size += block_size;
			rawData = reinterpret_cast<char*>(realloc(rawData, buf_size));
		}
		else {
			break;
		}
	} 

	rawData[rawSize] = 0; // add null termination

	status = true;
	return status;
}

 // *****************************************************************************************
 //
bool GenomeInputFile::load(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue,
	bool nonCanonical) {
	
	if (!status) {
		return false;
	}

	char* data;
	size_t total = 0;

	data = rawData;
	total = rawSize;
	
	if (status) {
		std::vector<char*> chromosomes;
		std::vector<size_t> lengths;
		std::vector<char*> headers;
		
		size_t totalLen = total;
		extractSubsequences(data, totalLen, chromosomes, lengths, headers);

		std::shared_ptr<MinHashFilter> minhashFilter = dynamic_pointer_cast<MinHashFilter>(filter);
		
		if (!minhashFilter) {
			throw std::runtime_error("unsupported filter type!");
		}
			
		kmerLength = minhashFilter->getLength();
		filterValue = minhashFilter->getFraction();
	
		// determine max k-mers count
		size_t sum_sizes = 0;
		for (auto e : lengths)
			sum_sizes += e - kmerLength + 1; 

		kmersBuffer.clear();
		kmersBuffer.resize(sum_sizes);

		kmersCount = 0;
		kmer_t* currentKmers = kmersBuffer.data();
		uint32_t* currentPositions = nullptr;

		if (storePositions) {
			positionsBuffer.clear();
			positionsBuffer.resize(sum_sizes);
			currentPositions = positionsBuffer.data();
		}

		for (size_t i = 0; i < chromosomes.size(); ++i) {
			size_t count = KmerHelper::extract(chromosomes[i], lengths[i], kmerLength, minhashFilter, currentKmers, currentPositions);
			currentKmers += count;
			currentPositions += storePositions ? count : 0;
			kmersCount += count;
		}
	
		kmers = kmersBuffer.data();

		//LOG_DEBUG << "Extraction: " << kmersCount << " kmers, " << chromosomes.size() << " chromosomes, " << totalLen << " bases" << endl ;
	}
	
	free(reinterpret_cast<void*>(rawData));
	
	return status;
}

// *****************************************************************************************
//
bool GenomeInputFile::initMultiFasta() {

	multifastaIndex = 0;

	if (!status) {
		return 0;
	}

	char* data;
	size_t total = 0;

	data = rawData;
	total = rawSize;

	size_t totalLen = total;
	extractSubsequences(data, totalLen, chromosomes, lengths, headers);
	
	return status;
}

// *****************************************************************************************
//
bool GenomeInputFile::loadNext(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue,
	bool nonCanonical,
	std::string& sampleName
) {
	
	std::shared_ptr<MinHashFilter> minhashFilter = dynamic_pointer_cast<MinHashFilter>(filter);

	if (!minhashFilter) {
		throw std::runtime_error("unsupported filter type!");
	}

	filterValue = minhashFilter->getFraction();
	kmerLength = minhashFilter->getLength();

	sampleName = headers[multifastaIndex];

	kmersBuffer.clear();
	kmersBuffer.resize(lengths[multifastaIndex]);

	kmersCount = KmerHelper::extract(
		chromosomes[multifastaIndex], 
		lengths[multifastaIndex], 
		kmerLength, 
		minhashFilter, 
		kmersBuffer.data(), 
		nullptr);

	kmers = kmersBuffer.data();
		
	++multifastaIndex;

	// no more sequences in multifasta
	return multifastaIndex < chromosomes.size();
}

// *****************************************************************************************
//
bool GenomeInputFile::extractSubsequences(
	char* data,
	size_t& totalLen,
	std::vector<char*>& subsequences,
	std::vector<size_t>& lengths,
	std::vector<char*>& headers) {

	// extract contigs
	char * header = nullptr;
	char * ptr = data;
	
	while ((header = strchr(ptr, '>'))) { // find begining of header
		*header = 0; // put 0 as separator (end of previous chromosome)
		if (subsequences.size()) {
			lengths.push_back(header - subsequences.back());
		}

		++header; // omit '<'
		headers.push_back(header);

		ptr = strchr(header, '\n'); // find end of header
		if (*(ptr - 1) == '\r') { // on Windows
			*(ptr - 1) = 0;
		}

		// !!! sd

		*ptr = 0; // put 0 as separator

		auto ptr_space = strchr(header, ' '); // find space in the header
		if (ptr_space)
			*ptr_space = 0;				// trim header to the first space

		++ptr; // move to next character (begin of chromosome)
		subsequences.push_back(ptr); // store chromosome
	}

	lengths.push_back(data + totalLen - subsequences.back());

	// remove newline characters from chromosome
	totalLen = 0;
	for (size_t i = 0; i < subsequences.size(); ++i) {
		// determine chromosome end
		char* newend = std::remove_if(subsequences[i], subsequences[i] + lengths[i], [](char c) -> bool { return c == '\n' || c == '\r';  });
		*newend = 0;
		lengths[i] = newend - subsequences[i];
		totalLen += lengths[i];
	//	assert(lengths[i] == strlen(subsequences[i]));
	}

	return true;
}

// *****************************************************************************************
//
bool MihashedInputFile::open(const std::string& filename)  {
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
	double& filterValue,
	bool nonCanonical) {
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

// *****************************************************************************************
//
bool KmcInputFile::load(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue,
	bool nonCanonical) {

	uint32_t counter;

	uint32 _mode;
	uint32 _counter_size;
	uint32 _lut_prefix_length;
	uint32 _signature_len;
	uint32 _min_count;
	uint64 _max_count;
	uint64 _total_kmers;

	kmcfile->Info(kmerLength, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	CKmerAPI kmer(kmerLength);

	uint64_t u_kmer;
	vector<uint64> tmp;

	// Wczytuje wszystkie k-mery z pliku do wektora, zeby pozniej moc robic prefetcha
	filter->setParams(kmerLength);
	std::shared_ptr<MinHashFilter> minhashFilter = std::dynamic_pointer_cast<MinHashFilter>(filter);

	// allocate buffers
#ifdef USE_RADULS
	kmersBuffer.resize(2 * _total_kmers + 4 * raduls::ALIGNMENT / sizeof(kmer_t));
	kmers = kmersBuffer.data();
	kmer_t* aux = kmersBuffer.data() + kmersBuffer.size() / 2;
	while (reinterpret_cast<uintptr_t>(kmers) % raduls::ALIGNMENT) ++kmers;
	while (reinterpret_cast<uintptr_t>(aux) % raduls::ALIGNMENT) ++aux;
#else
	kmersBuffer.resize(_total_kmers);
	kmers = kmersBuffer.data();
#endif

	// calculate k-mers shifting to get prefix of at least 8 bits
	size_t kmer_prefix_shift = 0;
	kmer_t tail_mask = 0;

	int prefix_bits = ((int)kmerLength - SUFFIX_LEN) * 2;

	if (prefix_bits < 8) {
		kmer_prefix_shift = (size_t)(8 - prefix_bits);
		tail_mask = (1ULL << kmer_prefix_shift) - 1;
	}

	kmersCount = 0;
	while (!kmcfile->Eof())
	{
		if (!kmcfile->ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();
		
		u_kmer = (u_kmer << kmer_prefix_shift) | (u_kmer & tail_mask);

		if ((*minhashFilter)(u_kmer)) {
			kmers[kmersCount++] = u_kmer;
		}
	}

#ifdef USE_RADULS
	size_t key_size = ((kmerLength * 2) + 7) / 8;
	raduls::PartialRadixSortMSD(reinterpret_cast<uint8_t*>(kmers), reinterpret_cast<uint8_t*>(aux), kmersCount, sizeof(kmer_t), key_size, key_size - 4, 4);
	if (key_size % 2) {
		std::swap(kmers, aux);
	}
#else
	//ParallelSort(kmers, kmersCount);
#endif

	filterValue = ((double)kmersCount / _total_kmers); // this may differ from theoretical
	//LOG_DEBUG << "Filter passed: " << kmersCount << "/" << _total_kmers << "(" << filterValue << ")" << endl ;
	filterValue = minhashFilter->getFraction(); // use proper value
	return kmcfile->Close();
}



