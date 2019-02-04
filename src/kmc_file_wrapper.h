#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>

class InputFile {
public:
	enum Format { KMC, MINHASH, GENOME };

	InputFile(std::shared_ptr<AbstractFilter> filter) : filter(filter) {}

	virtual bool open(const std::string& filename) = 0;
	virtual bool load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) = 0;
	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double filterValue) = 0;

protected:
	std::shared_ptr<AbstractFilter> filter;

};


class GenomeInputFile : public InputFile {
public:
	GenomeInputFile(std::shared_ptr<AbstractFilter> filter, bool storePositions)
		: InputFile(filter), status(false), isGzipped(false), storePositions(storePositions) {}

	virtual bool open(const std::string& filename) override;
	
	virtual bool load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) override;

	
	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double filterValue) override {
		throw std::runtime_error("Unable to store this kind of file");
	}

protected:
	size_t compressedSize;
	char* compressedData;

	bool status;
	bool isGzipped;
	bool storePositions;
	

};

class MihashedInputFile : public InputFile {
public:
	MihashedInputFile(std::shared_ptr<AbstractFilter> filter) : InputFile(filter), status(false) {}

	virtual bool open(const std::string& filename) override;
	
	virtual bool load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) override;

	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double filterValue) override;

protected:
	const uint32_t MINHASH_FORMAT_SIGNATURE = 0xfedcba98;

	std::vector<kmer_t> kmers;

	uint32_t kmerLength;

	double fraction;

	bool status;
};


class KmcInputFile : public InputFile {
public:
	
	KmcInputFile(std::shared_ptr<AbstractFilter> filter) : InputFile(filter) {}

	virtual bool open(const std::string& filename) override {
		kmcfile = std::make_shared<CKMCFile>();
		return kmcfile->OpenForListing(filename);
	}

	virtual bool load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) override;

	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double filterValue) override {
		throw std::runtime_error("Unable to store this kind of file");
	}

	
protected:
	std::shared_ptr<CKMCFile> kmcfile = nullptr;

};
