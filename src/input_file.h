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

	InputFile(std::shared_ptr<AbstractFilter> filter) : filter(filter) {}
	virtual ~InputFile() {};

	virtual bool open(const std::string& filename) = 0;
	virtual bool load(
		std::vector<kmer_t>& kmersBuffer, 
		std::vector<uint32_t>& positionsBuffer, 
		kmer_t*& kmers,
		size_t& kmersCount, 
		uint32_t& kmerLength, 
		double& filterValue,
		bool nonCanonical) = 0;


	static std::string removePathFromFile(const std::string& filePath) {
		size_t pos = filePath.find_last_of("/\\");
		if (pos != string::npos) {
			return filePath.substr(pos + 1);
		}
		else {
			return filePath;
		}
	}

protected:

	std::shared_ptr<AbstractFilter> filter;

};

// *****************************************************************************************
//
class GenomeInputFile : public InputFile {
public:
	GenomeInputFile(std::shared_ptr<AbstractFilter> filter, bool storePositions)
		: InputFile(filter), status(false), storePositions(storePositions) {}

	virtual ~GenomeInputFile() {}

	bool open(const std::string& filename) override;
	
	bool load(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue,
		bool nonCanonical) override;

	bool initMultiFasta();
	
	bool loadNext(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue,
		bool nonCanonical,
		std::string& sampleName,
		atomic<size_t>& total_kmers_in_kmers_collections
	);

	
protected:
	size_t rawSize;
	char* rawData;

	bool status;
	bool storePositions;

	// multifasta fields
	std::vector<char*> chromosomes;
	std::vector<size_t> lengths;
	std::vector<char*> headers;

	size_t multifastaIndex;
	
	bool extractSubsequences(
		char* data,
		size_t& totalLen,
		std::vector<char*>& subsequences,
		std::vector<size_t>& lengths,
		std::vector<char*>& headers);
};

// *****************************************************************************************
//
class MihashedInputFile : public InputFile {
public:
	MihashedInputFile(std::shared_ptr<AbstractFilter> filter) : InputFile(filter), status(false) {}
	virtual ~MihashedInputFile() {}

	bool open(const std::string& filename) override;
	
	bool load(
		std::vector<kmer_t>& kmersBuffer,
		std::vector<uint32_t>& positionsBuffer,
		kmer_t*& kmers,
		size_t& kmersCount,
		uint32_t& kmerLength,
		double& filterValue,
		bool nonCanonical) override;

	bool store(const std::string& filename, const kmer_t* kmers, size_t kmersCount, uint32_t kmerLength, double filterValue);

protected:
	const uint32_t MINHASH_FORMAT_SIGNATURE = 0xfedcba98;

	std::vector<kmer_t> kmers;

	uint32_t kmerLength;

	double fraction;

	bool status;
};

// *****************************************************************************************
//
class KmcInputFile : public InputFile {
public:
	
	KmcInputFile(std::shared_ptr<AbstractFilter> filter) : InputFile(filter) {}
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
		double& filterValue,
		bool nonCanonical) override;

protected:
	std::shared_ptr<CKMCFile> kmcfile{ nullptr };

};
