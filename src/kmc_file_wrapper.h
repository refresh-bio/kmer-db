#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"

#include <zlib.h>

#include <memory>
#include <fstream>
#include <string>
#include <cassert>


class InputFile {
public:
	enum Format { KMC, MINHASH, GENOME };

	InputFile(std::shared_ptr<IKmerFilter> filter) : filter(filter) {}

	virtual bool open(const std::string& filename) = 0;
	virtual bool load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& fraction) = 0;
	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) = 0;

protected:
	std::shared_ptr<IKmerFilter> filter;

};


class GenomeInputFile : public InputFile {
public:
	GenomeInputFile(std::shared_ptr<IKmerFilter> filter) : InputFile(filter), status(false) {}

	virtual bool open(const std::string& filename) override {
		
		status = false;
		gz = gzopen((filename + ".gz").c_str(), "rb"); // try .gz extension
		if (gz == Z_NULL) {
			gz = gzopen((filename + ".fasta").c_str(), "rb"); // try .fasta extension if .gz not found
		}

		if (gz != Z_NULL) {
			size_t blockSize = 1000000;
			data = reinterpret_cast<char*>(malloc(blockSize));
			char *ptr = data;

			total = 0;
			size_t got = 0;
			size_t allocated = blockSize;
			// read file in portions
			while ((got = gzread(gz, ptr, blockSize)) > 0) {
				total += got;
				allocated += blockSize;
				data = reinterpret_cast<char*>(realloc(data, allocated));
				ptr = data + total;
			}

			data[total] = 0;
			gzclose(gz);
			status = true;
		}

		return status;
	}

	virtual bool load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& fraction) override {
		if (!status) {
			return false;
		}

		std::vector<char*> chromosomes;
		std::vector<size_t> lengths;

		// extract contigs
		char * header = nullptr;
		char * ptr = data;

		while (header = strchr(ptr, '>')) { // find begining of header
			*header = 0; // put 0 as separator (end of previous chromosome)
			if (chromosomes.size()) {
				lengths.push_back(header - chromosomes.back());
			}

			++header;
			ptr = strchr(header, '\n'); // find end of header
			*ptr = 0; // put 0 as separator
			++ptr; // move to next character (begin of chromosome)
			chromosomes.push_back(ptr); // store chromosome
		}

		lengths.push_back(data + total - chromosomes.back());

		// remove newline characters from chromosome
		for (int i = 0; i < chromosomes.size(); ++i) {
			// determine chromosome end
			char* newend = std::remove_if(chromosomes[i], chromosomes[i] + lengths[i], [](char c) -> bool { return c == '\n' || c == '\r';  });
			*newend = 0;
			lengths[i] = newend - chromosomes[i];
			assert(lengths[i] == strlen(chromosomes[i]));
		}

		kmerLength = filter->getLength();

		// MAREK
		// Ekstrakcja k-merów z chromosomów. Wa¿ne zmienne:
		// - chromosomes - wektor ³añcuchów zakoñczonych zerem przechowuj¹cych sekwencje chromosomów
		// - lengths - wektor d³ugoœci chromosomów
		// - kmers - wektor wyjœciowy z kmerami
		// - filter - obiekt filtruj¹cy k-mery, operator () zwraca true jeœli k-mer spe³nia warunek
		// - kmerLength - d³ugoœæ k-mera

		

		// free memory
		free(reinterpret_cast<void*>(data));
	}

	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) override {
		throw std::runtime_error("Unable to store this kind of file");
	}

protected:
	char* data;
	gzFile gz;
	bool status;
	size_t total;

};

class MihashedInputFile : public InputFile {
public:
	MihashedInputFile(std::shared_ptr<IKmerFilter> filter) : InputFile(filter), status(false) {}

	virtual bool open(const std::string& filename) override {
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
	
	virtual bool load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& fraction) override {
		if (!status) {
			return false;
		}

		kmers = std::move(this->kmers);
		kmerLength = this->kmerLength;
		fraction = this->fraction;
		return true;
	}

	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) override {
		ofstream ofs(filename + ".minhash", std::ios_base::binary);
		ofs.write(reinterpret_cast<const char*>(&MINHASH_FORMAT_SIGNATURE), sizeof(MINHASH_FORMAT_SIGNATURE));
		size_t numKmers = kmers.size();
		ofs.write(reinterpret_cast<const char*>(&numKmers), sizeof(size_t));
		ofs.write(reinterpret_cast<const char*>(kmers.data()), kmers.size() * sizeof(kmer_t));
		ofs.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength));
		ofs.write(reinterpret_cast<const char*>(&fraction), sizeof(fraction));
		return true;
	}

protected:
	const uint32_t MINHASH_FORMAT_SIGNATURE = 0xfedcba98;

	std::vector<kmer_t> kmers;

	uint32_t kmerLength;

	double fraction;

	bool status;
};


class KmcInputFile : public InputFile {
public:
	
	KmcInputFile(std::shared_ptr<IKmerFilter> filter) : InputFile(filter) {}

	virtual bool open(const std::string& filename) override {
		kmcfile = std::make_shared<CKMCFile>();
		return kmcfile->OpenForListing(filename);
	}

	virtual bool load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& fraction) override {
		
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
		kmers.resize(_total_kmers);
		size_t passed = 0;
			
		filter->setParams(kmerLength);
			
		while (!kmcfile->Eof())
		{
			if (!kmcfile->ReadNextKmer(kmer, counter))
				break;
			kmer.to_long(tmp);
			u_kmer = tmp.front();

			if ((*filter)(u_kmer)) {
				kmers[passed++] = u_kmer;
			}
		}
		kmers.resize(passed);
		fraction = ((double)kmers.size() / _total_kmers); // this may differ from theoretical
		LOG_VERBOSE << "Min-hash passed: " << passed << "/" << _total_kmers << "(" << fraction << ")" << endl;
		return kmcfile->Close();
	}

	virtual bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double fraction) override {
		throw std::runtime_error("Unable to store this kind of file");
	}

	
protected:
	std::shared_ptr<CKMCFile> kmcfile = nullptr;

};
