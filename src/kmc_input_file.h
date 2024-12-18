#pragma once
#include "input_file.h"
#include "kmc_api/kmc_file.h"
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
template <class Filter>
class KmcInputFile : public FilteredInputFile<Filter> {
public:

	KmcInputFile(std::shared_ptr<Filter> filter) : FilteredInputFile<Filter>(filter) {}
	virtual ~KmcInputFile() {}

	bool open(const std::string& filename) override {
		kmcfile = std::make_shared<CKMCFile>();
		return kmcfile->OpenForListing(filename);
	}

	bool load(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue) override;

protected:
	std::shared_ptr<CKMCFile> kmcfile{ nullptr };

};


// *****************************************************************************************
//
template <class Filter>
bool KmcInputFile<Filter>::load(
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
	this->filter->configure(kmerLength);
	std::shared_ptr<MinHashFilter> minhashFilter = std::dynamic_pointer_cast<MinHashFilter>(this->filter);

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

	int prefix_bits = (int)kmerLength * 2 - SUFFIX_BITS;

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

		if ((*this->filter)(u_kmer)) {
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


