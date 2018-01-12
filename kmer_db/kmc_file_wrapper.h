#pragma once
#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"

#include <memory>
#include <fstream>

class KmcFileWrapper {
public:

	const uint32_t MINHASH_FORMAT_SIGNATURE = 0xfedcba98;

	KmcFileWrapper(std::shared_ptr<IKmerFilter> filter) : filter(filter) {}

	bool open(const std::string& filename, bool tryMinHash) {
		// try minhashed format
		if (tryMinHash) {
			auto file = std::make_shared<std::ifstream>(filename + ".minhash", std::ios_base::binary);
			if (*file) {
				uint32_t signature = 0;
				file->read(reinterpret_cast<char*>(&signature), sizeof(uint32_t));

				if (signature == MINHASH_FORMAT_SIGNATURE) {
					size_t numKmers;
					file->read(reinterpret_cast<char*>(&numKmers), sizeof(size_t));
					kmers.resize(numKmers);

					file->read(reinterpret_cast<char*>(kmers.data()), sizeof(kmer_t) * numKmers);
					file->read(reinterpret_cast<char*>(&kmerLength), sizeof(kmerLength));

					if (*file) {
						minhashFile = file;
						return true;
					}
					kmers.clear();
				}
			}
		}
		
		if (!minhashFile) {
			// try KMC format
			auto file = std::make_shared<CKMCFile>();
			if (file->OpenForListing(filename)) {
				kmcfile = file;
				return true;
			}
		}

		return false;
	}

	bool load(std::vector<kmer_t>& kmers, uint32_t& kmerLength) {
		if (minhashFile) {
			kmers = std::move(this->kmers);
			kmerLength = this->kmerLength;
			minhashFile->close();
			return true;
		}
		else if (kmcfile) {
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
			
			if (filter) {
				filter->setParams(kmerLength);
			}

			while (!kmcfile->Eof())
			{
				if (!kmcfile->ReadNextKmer(kmer, counter))
					break;
				kmer.to_long(tmp);
				u_kmer = tmp.front();

				if (filter == nullptr || (*filter)(u_kmer)) {
					kmers[passed++] = u_kmer;
				}
			}
			kmers.resize(passed);
			double fraction = (double)kmers.size() / _total_kmers;
			LOG_VERBOSE << "Min-hash passed: " << passed << "/" << _total_kmers << "(" << fraction << ")" << endl;
			return kmcfile->Close();
		}

		return false;
	}


	bool store(const std::string& filename, const std::vector<kmer_t>& kmers, uint32_t kmerLength) {
		ofstream ofs(filename + ".minhash", std::ios_base::binary);
		ofs.write(reinterpret_cast<const char*>(&MINHASH_FORMAT_SIGNATURE), sizeof(MINHASH_FORMAT_SIGNATURE));
		size_t numKmers = kmers.size();
		ofs.write(reinterpret_cast<const char*>(&numKmers), sizeof(size_t));
		ofs.write(reinterpret_cast<const char*>(kmers.data()), kmers.size() * sizeof(kmer_t));
		ofs.write(reinterpret_cast<const char*>(&kmerLength), sizeof(kmerLength));
		return true;
	}

protected:
	std::shared_ptr<CKMCFile> kmcfile = nullptr;
	std::shared_ptr<std::ifstream> minhashFile = nullptr;

	std::shared_ptr<IKmerFilter> filter;

	std::vector<kmer_t> kmers;

	uint32_t kmerLength;

};
