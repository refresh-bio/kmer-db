/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include "kmc_file_wrapper.h"
#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "kmer_extract.h"
#include "parallel_sorter.h"

#include <raduls.h>

#include <zlib.h>

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>


 bool GenomeInputFile::open(const std::string& filename) {

	status = false;

	FILE * in;
	
	if ((in = fopen((filename + ".gz").c_str(), "rb")) || 
		(in = fopen((filename + ".fna.gz").c_str(), "rb")) || 
		(in = fopen((filename + ".fasta.gz").c_str(), "rb"))) {
		isGzipped = true;
	}
	else {
		(in = fopen((filename + ".fna").c_str(), "rb")) || 
		(in = fopen((filename + ".fasta").c_str(), "rb"));
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

bool GenomeInputFile::load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) {
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
					ret = Z_STREAM_END;
					break;
				}

				if (ret == Z_OK && stream.avail_out == 0) {
					total = stream.total_out;
				}

				if (ret == Z_STREAM_END) {
					total = stream.total_out;
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
				data = reinterpret_cast<char*>(realloc(data, allocated));
				ptr = data + total;
			}

			inflateEnd(&stream);
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
		size_t totalLen = 0;

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
			totalLen += lengths[i];
			assert(lengths[i] == strlen(chromosomes[i]));
		}

		std::shared_ptr<MinHashFilter> minhashFilter = dynamic_pointer_cast<MinHashFilter>(filter);
		std::shared_ptr<SetFilter> setFilter = dynamic_pointer_cast<SetFilter>(filter);

		if (minhashFilter) {
			kmerLength = minhashFilter->getLength();
			filterValue = minhashFilter->getFilterValue();
		}
		else if (setFilter) {
		
		}
		else {
			throw std::runtime_error("unsupported filter type!");
		}

		//przewidywana liczba k-merow, zeby tylko raz byla alokacja pamieci w wektorze
		//przez 'N'ki i filtrowanie moze byc mniej faktycznie k-merow
		size_t sum_sizes = 0;
		for (auto e : lengths)
			sum_sizes += e - kmerLength + 1;

		kmers.clear();
		kmers.reserve(sum_sizes);

		if (storePositions) {
			positions.reserve(sum_sizes);
		}

		if (minhashFilter) {
			extractKmers(chromosomes, lengths, kmerLength, minhashFilter, kmers, positions, storePositions);
		}
		else {
			extractKmers(chromosomes, lengths, kmerLength, setFilter, kmers, positions, storePositions);
		}

		//auto bs = bitset<64>(kmers.front());

		//std::sort(kmers.begin(), kmers.end());
		ParallelSort(kmers.data(), kmers.size());
		auto it = std::unique(kmers.begin(), kmers.end());

		// iterate over kmers to select repeated ones
/*		size_t repeated = 0;
		auto nixt = std::next(kmers.begin());
		auto end = std::prev(kmers.end()); // the last kmer is either unique or will be reduced in std::unique operation
		for (auto it = kmers.begin(); it != end; ++it, ++nixt) {
			if (*it == *nixt) {
				*it |= KMER_MSB;
				++repeated;
			}
		}

	//	cout << "Repeated kmers: " << repeated << endl;
		auto it = std::unique(kmers.begin(), kmers.end(), 
			[](kmer_t a, kmer_t b)->bool { return (a << 1) == (b << 1);  }); // ignore MSB during comparison
	*/	
		kmers.erase(it, kmers.end());	

		LOG_DEBUG << "Extraction: " << kmers.size() << " kmers, " << chromosomes.size() << " chromosomes, " << totalLen << " bases" << endl;
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

bool MihashedInputFile::load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue) {
	if (!status) {
		return false;
	}

	kmers = std::move(this->kmers);
	kmerLength = this->kmerLength;
	filterValue = this->fraction;
	return true;
}

bool MihashedInputFile::store(const std::string& filename, const kmer_t* kmers, size_t kmersCount, uint32_t kmerLength, double filterValue) {
	ofstream ofs(filename + ".minhash", std::ios_base::binary);
	ofs.write(reinterpret_cast<const char*>(&MINHASH_FORMAT_SIGNATURE), sizeof(MINHASH_FORMAT_SIGNATURE));
	
	ofs.write(reinterpret_cast<const char*>(&kmersCount), sizeof(size_t));
	ofs.write(reinterpret_cast<const char*>(kmers), kmersCount * sizeof(kmer_t));
	ofs.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength));
	ofs.write(reinterpret_cast<const char*>(&filterValue), sizeof(filterValue));
	return true;
}




bool KmcInputFile::load(std::vector<kmer_t>& kmers, std::vector<uint32_t>& positions, uint32_t& kmerLength, double& filterValue)  {

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

	std::shared_ptr<MinHashFilter> minhashFilter = std::dynamic_pointer_cast<MinHashFilter>(filter);

	while (!kmcfile->Eof())
	{
		if (!kmcfile->ReadNextKmer(kmer, counter))
			break;
		kmer.to_long(tmp);
		u_kmer = tmp.front();

		if ((*minhashFilter)(u_kmer)) {
			kmers[kmersCount++] = u_kmer;
		}
	}
	kmers.resize(kmersCount * 2);

/*	auto prefix_comparer = [this](kmer_t a, kmer_t b)->bool {
		return GET_PREFIX(a) < GET_PREFIX(b);
	};
*/
	//std::sort(kmers.begin(), kmers.end(), prefix_comparer);
	ParallelSort(kmers.data(), kmers.size());


//	kmer_t* input = reinterpret_cast<uint8_t*>(kmers.data());
//	kmer_t* aux = kmers.data() + kmersCount;

//	raduls::RadixSortMSD(
//		 input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_threads);
	
	filterValue = ((double)kmers.size() / _total_kmers); // this may differ from theoretical
	LOG_DEBUG << "Filter passed: " << kmers.size() << "/" << _total_kmers << "(" << filterValue << ")" << endl;
	filterValue = minhashFilter->getFilterValue(); // use proper value
	return kmcfile->Close();
}



