#include "kmc_file_wrapper.h"

#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.1
Date   : 2018-06-12
*/

#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"

#include <zlib.h>

#include <memory>
#include <fstream>
#include <string>
#include <cassert>


 bool GenomeInputFile::open(const std::string& filename) {

	status = false;

	FILE * in;
	
	if ((in = fopen((filename + ".gz").c_str(), "rb")) || (in = fopen((filename + ".fna.gz").c_str(), "rb"))) {
		isGzipped = true;
	}
	else {
		in = fopen((filename + ".fna").c_str(), "rb");
	}
	
	if (in) {
		my_fseek(in, 0, SEEK_END);
		compressedSize = my_ftell(in);
		my_fseek(in, 0, SEEK_SET);

		compressedData = reinterpret_cast<char*>(malloc(compressedSize));
		fread(compressedData, compressedSize, 1, in);

		fclose(in);
		status = true;
	}

	return status;
}

bool GenomeInputFile::load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& filterValue) {
	if (!status) {
		return false;
	}

	char* data;
	
	size_t total = 0;

	if (isGzipped) {
		size_t blockSize = 10000000;
		data = reinterpret_cast<char*>(malloc(blockSize));

		// Init stream structure
		z_stream stream;
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
		stream.avail_in = compressedSize;
		stream.next_in = reinterpret_cast<Bytef*>(compressedData);

		if (inflateInit2(&stream, 31) == Z_OK) {

			// read data in portions
			char *ptr = data;

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
					status = false;
					ret == Z_STREAM_END;
					break;
				}

				if (ret == Z_STREAM_END) {
					break;
				}

				// reallocate only when some data left
				allocated += blockSize;
				data = reinterpret_cast<char*>(realloc(data, allocated));
				ptr = data + stream.total_out;
			}

			inflateEnd(&stream);
			total = stream.total_out;
		}
		else {
			status = false;
		}
	}
	else {
		data = compressedData;
		total = compressedSize;
	}

	if (status) {
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
		filterValue = filter->getFilterValue();

		// MAREK
		// Ekstrakcja k-merów z chromosomów. Wa¿ne zmienne:
		// - chromosomes - wektor ³añcuchów zakoñczonych zerem przechowuj¹cych sekwencje chromosomów
		// - lengths - wektor d³ugoœci chromosomów
		// - kmers - wektor wyjœciowy z kmerami
		// - filter - obiekt filtruj¹cy k-mery, operator () zwraca true jeœli k-mer spe³nia warunek
		// - kmerLength - d³ugoœæ k-mera

		static char* map = []() {
			static char _map[256];
			std::fill_n(_map, 256, -1);
			_map['a'] = _map['A'] = 0;
			_map['c'] = _map['C'] = 1;
			_map['g'] = _map['G'] = 2;
			_map['t'] = _map['T'] = 3;
			return _map;
		}();

		size_t total_kmers = 0;

		//przewidywana liczba k-merow, zeby tylko raz byla alokacja pamieci w wektorze
		//przez 'N'ki i filtrowanie moze byc mniej faktycznie k-merow
		size_t sum_sizes = 0;
		for (auto e : lengths)
			sum_sizes += e - kmerLength + 1;

		//filter->initialize(sum_sizes);
		size_t kmersCount = 0;
		kmers.clear();
		kmers.reserve(sum_sizes);

		kmer_t kmer_str, kmer_rev, kmer_can;
		uint32_t kmer_len_shift = (kmerLength - 1) * 2;
		kmer_t kmer_mask = (1ull << (2 * kmerLength)) - 1;
		int omit_next_n_kmers;
		uint32_t i;
		for (size_t j = 0; j < chromosomes.size(); ++j)
		{
			char* seq = chromosomes[j];
			size_t seq_size = lengths[j];

			kmer_str = kmer_rev = 0;

			uint32_t str_pos = kmer_len_shift - 2;
			uint32_t rev_pos = 2;

			omit_next_n_kmers = 0;

			for (i = 0; i < kmerLength - 1; ++i, str_pos -= 2, rev_pos += 2)
			{
				char symb = map[seq[i]];
				if (symb < 0)
				{
					symb = 0;
					omit_next_n_kmers = i + 1;
				}
				kmer_str += (kmer_t)symb << str_pos;
				kmer_rev += (kmer_t)(3 - symb) << rev_pos;
			}

			for (; i < seq_size; ++i)
			{
				char symb = map[seq[i]];
				if (symb < 0)
				{
					symb = 0;
					omit_next_n_kmers = kmerLength;
				}
				kmer_str = (kmer_str << 2) + (kmer_t)symb;
				kmer_str &= kmer_mask;

				kmer_rev >>= 2;
				kmer_rev += (kmer_t)(3 - symb) << kmer_len_shift;

				if (omit_next_n_kmers > 0)
				{
					--omit_next_n_kmers;
					continue;
				}

				kmer_can = (kmer_str < kmer_rev) ? kmer_str : kmer_rev;

				//filter->add(kmer_can);	
				if ((*filter)(kmer_can)) {
					kmers.push_back(kmer_can);
				}
			}
		}
		//kmers = std::move(filter->finalize());

		std::sort(kmers.begin(), kmers.end());
		auto it = std::unique(kmers.begin(), kmers.end());
		kmers.erase(it, kmers.end());
	}
	
	// free memory
	if (data != compressedData) {
		free(reinterpret_cast<void*>(data));
	}
	free(reinterpret_cast<void*>(compressedData));
	
	return status;
}


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

bool MihashedInputFile::load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& filterValue) {
	if (!status) {
		return false;
	}

	kmers = std::move(this->kmers);
	kmerLength = this->kmerLength;
	filterValue = this->fraction;
	return true;
}

bool MihashedInputFile::store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength, double filterValue) {
	ofstream ofs(filename + ".minhash", std::ios_base::binary);
	ofs.write(reinterpret_cast<const char*>(&MINHASH_FORMAT_SIGNATURE), sizeof(MINHASH_FORMAT_SIGNATURE));
	size_t numKmers = kmers.size();
	ofs.write(reinterpret_cast<const char*>(&numKmers), sizeof(size_t));
	ofs.write(reinterpret_cast<const char*>(kmers.data()), kmers.size() * sizeof(kmer_t));
	ofs.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength));
	ofs.write(reinterpret_cast<const char*>(&filterValue), sizeof(filterValue));
	return true;
}




bool KmcInputFile::load(std::vector<kmer_t>& kmers, uint32_t& kmerLength, double& filterValue)  {

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
	kmers.resize(_total_kmers);
	size_t kmersCount = 0;
	//filter->initialize(_total_kmers);

	size_t passed = 0;

	while (!kmcfile->Eof())
	{
		if (!kmcfile->ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		if ((*filter)(u_kmer)) {
			kmers[kmersCount++] = u_kmer;
		}
	}
	kmers.resize(kmersCount);
	//kmers = std::move(filter->finalize());

	filterValue = ((double)kmers.size() / _total_kmers); // this may differ from theoretical
	LOG_VERBOSE << "Filter passed: " << kmers.size() << "/" << _total_kmers << "(" << filterValue << ")" << endl;
	return kmcfile->Close();
}



