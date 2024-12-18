#pragma once
#include "input_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "kmer_extract.h"
#include "alphabet.h"

#include "../libs/refresh/compression/lib/file_wrapper.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>
#include <filesystem>
#include <memory>


// *****************************************************************************************
//
template <class Filter>
class GenomeInputFile : public FilteredInputFile<Filter>, public IMultiSampleFile {
public:
	GenomeInputFile(std::shared_ptr<Filter> filter, std::shared_ptr<Alphabet> alphabet)
		: FilteredInputFile<Filter>(filter), alphabet(alphabet) {}

	virtual ~GenomeInputFile() {}

	bool open(const std::string& filename) override;

	bool load(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue) override;

	bool initMultiFasta() override;

	bool loadNext(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue,
		std::string& sampleName,
		atomic<size_t>& total_kmers_in_kmers_collections) override;


protected:
	size_t rawSize{ 0 };
	char* rawData{ nullptr };

	bool status{ false };
	bool storePositions{ false };

	size_t multifastaIndex{ 0 };

	// multifasta fields
	std::vector<char*> chromosomes;
	std::vector<size_t> lengths;
	std::vector<char*> headers;

	std::shared_ptr<Alphabet> alphabet;

	inline bool extractSubsequences(
		char* data,
		size_t& totalLen,
		std::vector<char*>& subsequences,
		std::vector<size_t>& lengths,
		std::vector<char*>& headers);
};


// *****************************************************************************************
//
template <class Filter>
bool GenomeInputFile<Filter>::open(const std::string& filename) {

	vector<string> extensions{
		"", ".fa", ".fna", ".fasta", 
		".gz", ".fa.gz", ".fna.gz", ".fasta.gz" };

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
template <class Filter>
bool GenomeInputFile<Filter>::load(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue) {

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
		
		std::shared_ptr<MinHashFilter> minhashFilter = std::dynamic_pointer_cast<MinHashFilter>(this->filter);

		if (!minhashFilter) {
			throw std::runtime_error("Only MinHashFilter is currently supported!");
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

			size_t count = KmerHelper::extract<Filter>(
				chromosomes[i], lengths[i], kmerLength, *this->alphabet, *this->filter, currentKmers);
			
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
template <class Filter>
bool GenomeInputFile<Filter>::initMultiFasta() {

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
template <class Filter>
bool GenomeInputFile<Filter>::loadNext(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue,
	std::string& sampleName,
	atomic<size_t>& total_kmers_in_kmers_collections
) {

	std::shared_ptr<MinHashFilter> minhashFilter = dynamic_pointer_cast<MinHashFilter>(this->filter);

	if (!minhashFilter) {
		throw std::runtime_error("Only MinHashFilter is currently supported!");
	}

	filterValue = minhashFilter->getFraction();
	kmerLength = minhashFilter->getLength();

	sampleName = headers[multifastaIndex];

	kmersBuffer.clear();
	total_kmers_in_kmers_collections -= kmersBuffer.size();
	kmersBuffer.resize(lengths[multifastaIndex]);
	total_kmers_in_kmers_collections += kmersBuffer.size();

	kmersCount = KmerHelper::extract<Filter>(
		chromosomes[multifastaIndex],
		lengths[multifastaIndex],
		kmerLength,
		*this->alphabet,
		*this->filter,
		kmersBuffer.data());
	
	kmers = kmersBuffer.data();

	++multifastaIndex;

	// no more sequences in multifasta
	return multifastaIndex < chromosomes.size();
}

// *****************************************************************************************
//
template <class Filter>
bool GenomeInputFile<Filter>::extractSubsequences(
	char* data,
	size_t& totalLen,
	std::vector<char*>& subsequences,
	std::vector<size_t>& lengths,
	std::vector<char*>& headers) {

	// extract contigs
	char* header = nullptr;
	char* ptr = data;

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