#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "loader_tasks.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>

// *****************************************************************************************
//
class InputFile {
public:
	enum Format { KMC, MINHASH, GENOME };

	static std::string format2string(enum Format f) {
		switch (f) {
		case GENOME: return "fasta genomes";
		case KMC: return "k-mers";
		case MINHASH: return "minhashed k-mers";
		}

		return "";
	}

	virtual bool open(const std::string& filename) = 0;
	
	virtual bool load(
		std::vector<kmer_t>& kmersBuffer, 
		std::vector<uint32_t>& positionsBuffer, 
		kmer_t*& kmers,
		size_t& kmersCount, 
		uint32_t& kmerLength, 
		double& filterValue) = 0;

	virtual ~InputFile() {}
};


//******************************************************************************************
//
class IMultiSampleFile {

public:

	virtual bool initMultiFasta() = 0;

	virtual bool loadNext(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue,
		std::string& sampleName,
		atomic<size_t>& total_kmers_in_kmers_collections) = 0;

	virtual ~IMultiSampleFile() {}
};


// *****************************************************************************************
//
template <class Filter>
class FilteredInputFile : public InputFile {
public:
	FilteredInputFile(std::shared_ptr<Filter>& filter) : filter(filter) {}

protected:
	std::shared_ptr<Filter> filter;
};


