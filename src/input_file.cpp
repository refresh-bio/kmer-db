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

#include <zlib.h>

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>

#ifdef USE_RADULS
	#include <raduls.h>
#endif

// *****************************************************************************************
//
 bool GenomeInputFile::open(const std::string& filename) {

	status = false;

	FILE * in;
	
	// try to open without adding extension
	if ((in = fopen(filename.c_str(), "rb"))) {
		isGzipped = filename.substr(filename.length() - 3) == ".gz";
	}
	else {
		// try adding an extension
		if (
			(in = fopen((filename + ".gz").c_str(), "rb")) ||
			(in = fopen((filename + ".fa.gz").c_str(), "rb")) ||
			(in = fopen((filename + ".fna.gz").c_str(), "rb")) ||
			(in = fopen((filename + ".fasta.gz").c_str(), "rb"))) {
			isGzipped = true;
		}
		else {
			(in = fopen((filename + ".fa").c_str(), "rb")) ||
			(in = fopen((filename + ".fna").c_str(), "rb")) ||
			(in = fopen((filename + ".fasta").c_str(), "rb"));
		}
	}

	if (in) {
		my_fseek(in, 0, SEEK_END);
		rawSize = my_ftell(in);
		my_fseek(in, 0, SEEK_SET);

		rawData = reinterpret_cast<char*>(malloc(rawSize + 1));
		size_t blocksRead = fread(rawData, rawSize, 1, in);
		rawData[rawSize] = 0; // add null termination 
		fclose(in);
		status = blocksRead;
	}

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
	double& filterValue) {
	
	if (!status) {
		return false;
	}

	char* data;
	size_t total = 0;

	if (isGzipped) {
		status = unzip(rawData, rawSize, data, total);
	}
	else {
		data = rawData;
		total = rawSize;
	}

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
			size_t count = extractKmers(chromosomes[i], lengths[i], kmerLength, minhashFilter, currentKmers, currentPositions);
			currentKmers += count;
			currentPositions += storePositions ? count : 0;
			kmersCount += count;
		}
	
		ParallelSort(kmersBuffer.data(), kmersCount);
		auto it = std::unique(kmersBuffer.begin(), kmersBuffer.begin() + kmersCount);
		
		kmers = kmersBuffer.data();
		kmersCount = it - kmersBuffer.begin();

		LOG_DEBUG << "Extraction: " << kmersCount << " kmers, " << chromosomes.size() << " chromosomes, " << totalLen << " bases" << endl;
	}
	
	// free memory
	if (data != rawData) {
		free(reinterpret_cast<void*>(data));
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

	if (isGzipped) {
		status = unzip(rawData, rawSize, data, total);
	}
	else {
		data = rawData;
		total = rawSize;
	}

	if (status) {
		size_t totalLen = total;
		extractSubsequences(data, totalLen, chromosomes, lengths, headers);
	}

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
	std::string& sampleName
) {
	
	// no more sequences in multifasta
	if (multifastaIndex == chromosomes.size()) {
		return false;
	}

	std::shared_ptr<MinHashFilter> minhashFilter = dynamic_pointer_cast<MinHashFilter>(filter);

	if (!minhashFilter) {
		throw std::runtime_error("unsupported filter type!");
	}

	filterValue = minhashFilter->getFraction();
	kmerLength = minhashFilter->getLength();
	sampleName = headers[multifastaIndex];

	kmersBuffer.clear();
	kmersBuffer.resize(lengths[multifastaIndex]);

	kmersCount = extractKmers(
		chromosomes[multifastaIndex], 
		lengths[multifastaIndex], 
		kmerLength, 
		minhashFilter, 
		kmersBuffer.data(), 
		nullptr);

	ParallelSort(kmersBuffer.data(), kmersCount);
	auto it = std::unique(kmersBuffer.begin(), kmersBuffer.begin() + kmersCount);

	kmers = kmersBuffer.data();
	kmersCount = it - kmersBuffer.begin();
		
	++multifastaIndex;

	return true;
}


// *****************************************************************************************
//
bool GenomeInputFile::unzip(char* compressedData, size_t compressedSize, char*&outData, size_t &outSize) {
	bool ok = true;
	size_t blockSize = 10000000;
	outData = reinterpret_cast<char*>(malloc(blockSize));

	// Init stream structure
	z_stream stream;
	stream.zalloc = Z_NULL;
	stream.zfree = Z_NULL;
	stream.opaque = Z_NULL;
	stream.avail_in = compressedSize;
	stream.next_in = reinterpret_cast<Bytef*>(compressedData);

	if (inflateInit2(&stream, 31) == Z_OK) {

		// read data in portions
		char *ptr = outData;

		size_t allocated = blockSize;

		// decompress file in portions
		for (;;) {
			stream.avail_out = blockSize;
			stream.next_out = reinterpret_cast<Bytef*>(ptr);
			int ret = inflate(&stream, Z_NO_FLUSH);

			switch (ret)
			{
			case Z_NEED_DICT:
			case Z_DATA_ERROR:
			case Z_MEM_ERROR:
				ok = false;
				ret = Z_STREAM_END;
				break;
			}

			if (ret == Z_OK && stream.avail_out == 0) {
				outSize = stream.total_out;
			}

			if (ret == Z_STREAM_END) {
				outSize = stream.total_out;
				//multistream detection
				if (stream.avail_in >= 2 && stream.next_in[0] == 0x1f && stream.next_in[1] == 0x8b) {
					if (inflateReset(&stream) != Z_OK) {
						LOG_NORMAL << "Error while reading gzip file\n";
						exit(1);
					}
				}
				else
					break;
			}

			// reallocate only when some data left
			allocated += blockSize;
			outData = reinterpret_cast<char*>(realloc(outData, allocated + 1)); // allocate for null termination
			ptr = outData + outSize;
		}

		inflateEnd(&stream);
		outData[outSize] = 0;
	}
	else {
		ok = false;
	}

	return ok;
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
		*ptr = 0; // put 0 as separator
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

// *****************************************************************************************
//
bool KmcInputFile::load(
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer,
	kmer_t*& kmers,
	size_t& kmersCount,
	uint32_t& kmerLength,
	double& filterValue) {

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
	ParallelSort(kmers, kmersCount);
#endif

	filterValue = ((double)kmersCount / _total_kmers); // this may differ from theoretical
	LOG_DEBUG << "Filter passed: " << kmersCount << "/" << _total_kmers << "(" << filterValue << ")" << endl;
	filterValue = minhashFilter->getFraction(); // use proper value
	return kmcfile->Close();
}



