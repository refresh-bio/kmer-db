#include "params.h"

#include <iostream>
#include <stdexcept>

using namespace std;


// *****************************************************************************************
//
Params::Params() {
	availableMetrics["jaccard"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / (cnt1 + cnt2 - common); };
	availableMetrics["min"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::min(cnt1, cnt2); };
	availableMetrics["max"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / std::max(cnt1, cnt2); };
	availableMetrics["cosine"] = [](size_t common, size_t cnt1, size_t cnt2, int kmerLength) -> double { return (double)common / sqrt(cnt1 * cnt2); };
	availableMetrics["mash"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / (queryCnt + dbCnt - common);
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
	};
	availableMetrics["ani"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / (queryCnt + dbCnt - common);
		double d_mash = (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
		return 1.0 - d_mash;
	};
	availableMetrics["ani-shorter"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / min(queryCnt, dbCnt);
		double d_mash = (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
		return 1.0 - d_mash;
	};

	availableMetrics["mash-query"] = [](size_t common, size_t queryCnt, size_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / queryCnt;
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
	};
}

// *****************************************************************************************
//
Params::Mode Params::str2mode(const std::string& str) {
	if (str == MODE_BUILD) { return Mode::build; }
	else if (str == MODE_MINHASH) { return Mode::minhash; }
	else if (str == MODE_ALL_2_ALL) { return Mode::all2all; }
	else if (str == MODE_ALL_2_ALL_SPARSE) { return Mode::all2all_sparse; }
	else if (str == MODE_ALL_2_ALL_PARTS) { return Mode::all2all_parts; }
	else if (str == MODE_NEW_2_ALL) { return Mode::new2all; }
	else if (str == MODE_ONE_2_ALL) { return Mode::one2all; }
	else if (str == MODE_DISTANCE) { return Mode::distance; }
	else { return Mode::unknown; }
}

// *****************************************************************************************
//
bool Params::parse(int argc, char** argv) {

	std::vector<string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* c)->string { return c; });

	bool helpWanted = findSwitch(params, SWITCH_HELP);
	if (params.size() == 0) {
		// no arguments or -help switch already consumed
		showInstructions(Mode::unknown);
		return false;
	}
	else if (helpWanted && params.size() == 1) {
		// help for particular mode
		showInstructions(str2mode(params[0]));
		return false;
	}
	else {

		// search for switches and options
		if (findSwitch(params, OPTION_VERBOSE)) { // verbose mode
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
		}

		if (findSwitch(params, OPTION_DEBUG)) { // verbose mode
			Log::getInstance(Log::LEVEL_VERBOSE).enable();
			Log::getInstance(Log::LEVEL_DEBUG).enable();
		}

		multisampleFasta = findSwitch(params, SWITCH_MULTISAMPLE_FASTA);
		
		bool fractionSpecified = findOption(params, OPTION_FRACTION, fraction);				// minhash fraction
		findOption(params, OPTION_FRACTION_START, fractionStart);	// minhash fraction start value
		findOption(params, OPTION_LENGTH, kmerLength);				// kmer length

		if (kmerLength > 30) {
			throw std::runtime_error("K-mer length cannot not exceed 30.");
		}

		findOption(params, OPTION_THREADS, numThreads);			// number of threads
		if (numThreads <= 0) {
			numThreads = std::max((int)std::thread::hardware_concurrency(), 1); // hardware_concurrency may return 0
		}

		findOption(params, OPTION_READER_THREADS, numReaderThreads);	// number of threads
		if (numReaderThreads <= 0) {
			// more reader threads for smaller filters (from t/8 up to t)
			int invFraction = (int)(1.0 / fraction);
			numReaderThreads = std::max(std::min(numThreads, (numThreads / 8) * invFraction), 1);
		}

		findOption(params, OPTION_BUFFER, cacheBufferMb);	// size of temporary buffer in megabytes
		if (cacheBufferMb <= 0) {
			cacheBufferMb = 8;
		}

		findOption(params, OPTION_BELOW, below);
		findOption(params, OPTION_ABOVE, above);

		bool below_eq = findOption(params, OPTION_BELOW_EQ, below);
		bool above_eq = findOption(params, OPTION_ABOVE_EQ, above);

		if (findSwitch(params, SWITCH_KMC_SAMPLES)) {
			inputFormat = InputFile::KMC;
			kmerLength = 0;
		}

		if (findSwitch(params, SWITCH_MINHASH_SAMPLES)) {
			if (inputFormat == InputFile::KMC) {
				throw std::runtime_error(
					SWITCH_KMC_SAMPLES + " and " + SWITCH_MINHASH_SAMPLES + " switches exclude one another.");
			}
			inputFormat = InputFile::MINHASH;
			fraction = 1.0;
			kmerLength = 0;
		}

		sparseOut = findSwitch(params, SWITCH_SPARSE);
		extendDb = findSwitch(params, SWITCH_EXTEND_DB);
		nonCanonical = findSwitch(params, SWITCH_NON_CANONICAL);

		phylipOut = findSwitch(params, SWITCH_PHYLIP_OUT);
		if (phylipOut) {
			sparseOut = false;
		}

		if (params.size()) {
			string mode = params.front();
			params.erase(params.begin());

			// detect obsolete modes
			if (mode == "build-kmers" || mode == "build-mh") {
				throw std::runtime_error(
					"build -kmers / build -mh modes are obsolete, use " + SWITCH_KMC_SAMPLES + " / " + SWITCH_MINHASH_SAMPLES + " switches instead.");
			}

			this->mode = str2mode(mode);

			// fill metrics in distance mode
			if (this->mode == Mode::distance) {
				for (const auto& entry : availableMetrics) {
					if (findSwitch(params, entry.first)) {
						metrics.push_back(entry);
					}
				}

				// if empty, add jacard
				if (metrics.empty()) {
					metrics.push_back(*availableMetrics.find("jaccard"));
				}

				// alter thresholds if above_eq or below_eq are selected 
				if (below_eq) {
					below = std::nextafter(below, below + 1.0);
				}

				if (above_eq) {
					above = std::nextafter(above, above - 1.0);
				}
			}
			else {
				// alter thresholds if above_eq or below_eq are selected 
				if (below_eq) {
					below += 1.0;
				}

				if (above_eq) {
					above -= 1.0;
				}
			}

			// set default fration in minhash mode if not specified
			if (this->mode == Mode::minhash && !fractionSpecified) {
				fraction = 0.01;
			}

			files.insert(files.begin(), params.begin(), params.end());
		}
		else {
			throw std::runtime_error("Mode not specified");
		}
	}

	return true;
}

// *****************************************************************************************
//
void Params::showInstructions(Mode mode) const {

	if (mode == Mode::build) {
		LOG_NORMAL
			<< "Building a database:" << endl
			<< "    kmer-db " << MODE_BUILD
			<< " [" << OPTION_LENGTH << " <kmer-length>]"
			<< " [" << OPTION_FRACTION << " <fraction>]"
			<< " [" << SWITCH_MULTISAMPLE_FASTA << "]"
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_THREADS << " <threads>] <sample_list> <database>" << endl

			<< "    kmer-db " << MODE_BUILD << " " << SWITCH_KMC_SAMPLES
			<< " [" << OPTION_FRACTION << " <fraction>]"
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_THREADS << " <threads>] <sample_list> <database>" << endl

			<< "    kmer-db " << MODE_BUILD << " " << SWITCH_MINHASH_SAMPLES
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_THREADS << " <threads>] <sample_list> <database>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    sample_list (input) - file containing list of samples in one of the following formats:" << endl
			<< "                          FASTA genomes/reads (default), KMC k-mers (" << SWITCH_KMC_SAMPLES << "), or minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << ")," << endl
			<< "    database (output) - file with generated k-mer database," << endl

			<< "Options: " << endl
			<< "    " << OPTION_LENGTH << " <kmer_length> - length of k-mers (default: 18, maximum: 30)" << endl
			<< "    " << OPTION_FRACTION << " <fraction> - fraction of all k-mers to be accepted by the minhash filter (default: 1)" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl
			<< "    " << SWITCH_EXTEND_DB << " - extend the existing database with new samples" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl << endl;
	}
	else if (mode == Mode::all2all) {
		LOG_NORMAL
			<< "Counting common k-mers for all the samples in the database:" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << SWITCH_SPARSE << " [" << OPTION_ABOVE << " <v>] [" << OPTION_BELOW << " <v>] [" << OPTION_ABOVE_EQ << " <v>] [" << OPTION_BELOW_EQ << " <v>]]"
			<< " <database> <common_table>" << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << SWITCH_SPARSE << " - produce sparse matrix as output" << endl
			<< "    " << OPTION_ABOVE << " <v> - retains elements larger then <v>" << endl
			<< "    " << OPTION_BELOW << " <v> - retains elements smaller then <v>" << endl
			<< "    " << OPTION_ABOVE_EQ << " <v> - retains elements larger or equal <v>" << endl
			<< "    " << OPTION_BELOW_EQ << " <v> - retains elements smaller or equal <v>" << endl
			<< endl;
	}
	else if (mode == Mode::all2all_sparse) {
		LOG_NORMAL
			<< "Counting common k-mers for all the samples in the database (sparse computation):" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL_SPARSE
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << OPTION_ABOVE << " <v>] [" << OPTION_BELOW << " <v>] [" << OPTION_ABOVE_EQ << " <v>] [" << OPTION_BELOW_EQ << " <v>]"
			<< " <database> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << OPTION_ABOVE << " <v> - retains elements larger then <v>" << endl
			<< "    " << OPTION_BELOW << " <v> - retains elements smaller then <v>" << endl
			<< "    " << OPTION_ABOVE_EQ << " <v> - retains elements larger or equal <v>" << endl
			<< "    " << OPTION_BELOW_EQ << " <v> - retains elements smaller or equal <v>" << endl
			<< endl;
	}
	else if (mode == Mode::all2all_parts) {
		LOG_NORMAL
			<< "Counting common k-mers for all the samples in the database parts (sparse computation):" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL_PARTS
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << OPTION_ABOVE << " <v>] [" << OPTION_BELOW << " <v>] [" << OPTION_ABOVE_EQ << " <v>] [" << OPTION_BELOW_EQ << " <v>]"
			<< " <db_list> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    db_list (input) - file containing list of database file names" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << OPTION_ABOVE << " <v> - retains elements larger then <v>" << endl
			<< "    " << OPTION_BELOW << " <v> - retains elements smaller then <v>" << endl
			<< "    " << OPTION_ABOVE_EQ << " <v> - retains elements larger or equal <v>" << endl
			<< "    " << OPTION_BELOW_EQ << " <v> - retains elements smaller or equal <v>" << endl
			<< endl;
	}
	/*
	else if (mode == Mode::db2db) {
		 LOG_NORMAL
			 << "Counting common k-mers for all the samples in between the databases (sparse computation):" << endl
			 << "    kmer-db " << MODE_DB_2_DB_SP
			 << " [" << OPTION_BUFFER << " <size_mb>]"
			 << " [" << SWITCH_SPARSE << "]"
			 << " [" << OPTION_THREADS << " <threads>] <database> <common_table>" << endl << endl

			 << "Positional arguments:" << endl
			 << "    database1 (input) - k-mer database file" << endl
			 << "    database2 (input) - k-mer database file" << endl
			 << "    common_table (output) - comma-separated table with number of common k-mers" << endl

			 << "Options:" << endl
			 << "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			 << "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			 << "    " << SWITCH_SPARSE << " - produce sparse matrix as output" << endl
			 << "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl << endl;
	 }
	 */
	else if (mode == Mode::new2all) {
		LOG_NORMAL
			<< "Counting common kmers between set of new samples and all the samples in the database:" << endl
			<< "    kmer-db " << MODE_NEW_2_ALL
			<< " [" << SWITCH_MULTISAMPLE_FASTA << " | " << SWITCH_KMC_SAMPLES << " | " << SWITCH_MINHASH_SAMPLES << "]"
			<< " [" << SWITCH_SPARSE << " [" << OPTION_ABOVE << " <v>] [" << OPTION_BELOW << " <v>] [" << OPTION_ABOVE_EQ << " <v>] [" << OPTION_BELOW_EQ << " <v>]] "
			<< " [" << OPTION_THREADS << " <threads>] <database> <sample_list> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    sample_list (input) - file containing list of samples in one of the following formats:" << endl
			<< "                          FASTA genomes/reads (default), KMC k-mers (" << SWITCH_KMC_SAMPLES << "), or minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << ")" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << SWITCH_SPARSE << " - outputs a sparse matrix" << endl
			<< "    " << OPTION_ABOVE << " <v> - retains elements larger then <v>" << endl
			<< "    " << OPTION_BELOW << " <v> - retains elements smaller then <v>" << endl
			<< "    " << OPTION_ABOVE_EQ << " <v> - retains elements larger or equal <v>" << endl
			<< "    " << OPTION_BELOW_EQ << " <v> - retains elements smaller or equal <v>" << endl
			<< endl;
	}
	else if (mode == Mode::one2all) {
		LOG_NORMAL
			<< "Counting common kmers between single sample and all the samples in the database:" << endl
			<< "    kmer-db " << MODE_ONE_2_ALL
			<< " [" << SWITCH_KMC_SAMPLES << " | " << SWITCH_MINHASH_SAMPLES << "] <database> <sample> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file." << endl
			<< "    sample (input) - query sample in one of the supported formats:" << endl
			<< "                     FASTA genomes/reads (default), KMC k-mers (" << SWITCH_KMC_SAMPLES << "), or minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << "), " << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers." << endl << endl;
	}
	else if (mode == Mode::distance) {
		LOG_NORMAL
			<< "Calculating similarities/distances on the basis of common k-mers:" << endl
			<< "    kmer-db " << MODE_DISTANCE << " [<measures>]"
			<< " [" << SWITCH_SPARSE << " [" << OPTION_ABOVE << " <v>] [" << OPTION_BELOW << " <v>] [" << OPTION_ABOVE_EQ << " <v>] [" << OPTION_BELOW_EQ << " <v>]] "
			<< "<common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    common_table (input) - comma-separated table with a number of common k-mers" << endl

			<< "Options:" << endl
			<< "    measures - names of the similarity/distance measures to be calculated, one or more of the following" << endl
			<< "               jaccard (default), min, max, cosine, mash, ani, ani-shorter." << endl
			<< "    " << SWITCH_PHYLIP_OUT << " - store output distance matrix in a Phylip format" << endl
			<< "    " << SWITCH_SPARSE << " - outputs a sparse matrix (independently of the input matrix format)" << endl
			<< "    " << OPTION_ABOVE << " <v> - retains elements larger then <v>" << endl
			<< "    " << OPTION_BELOW << " <v> - retains elements smaller then <v>" << endl 
			<< "    " << OPTION_ABOVE_EQ << " <v> - retains elements larger or equal <v>" << endl
			<< "    " << OPTION_BELOW_EQ << " <v> - retains elements smaller or equal <v>" << endl
			<< endl
			<< "This mode generates a file with similarity/distance table for each selected measure." << endl
			<< "Name of the output file is produced by adding to the input file an extension with a measure name." << endl << endl;
	}
	else if (mode == Mode::minhash) {
		LOG_NORMAL
			<< "Storing minhashed k-mers:" << endl
			<< "    kmer-db " << MODE_MINHASH << " [" << OPTION_FRACTION << " <fraction> ][" << OPTION_LENGTH << " <kmer_length>][" << SWITCH_MULTISAMPLE_FASTA << "] <sample_list>" << endl
			<< "    kmer-db " << MODE_MINHASH << " " << SWITCH_KMC_SAMPLES << " [" << OPTION_FRACTION << " <fraction>] <sample_list>" << endl << endl
			<< "Positional arguments:" << endl
			<< "    sample (input) - query sample in one of the supported formats:" << endl
			<< "                     FASTA genomes/reads (default) or KMC k-mers (" << SWITCH_KMC_SAMPLES << ")" << endl
			<< "Options:" << endl
			<< "    " << OPTION_FRACTION << " <fraction> - fraction of all k-mers to be accepted by the minhash filter (default: 0.01)" << endl
			<< "    " << OPTION_LENGTH << " <kmer_length> - length of k-mers (default: 18, maximum: 30)" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl << endl
			<< "For each sample from the list, a binary file with *.minhash* extension containing filtered k-mers is created" << endl << endl;
	}
	else {
		LOG_NORMAL
			<< "USAGE" << endl
			<< "    kmer-db <mode> [options] <positional arguments>" << endl << endl

			<< "Modes:" << endl
			<< "    " << MODE_BUILD << " - building a database from FASTA genomes/reads, k-mers, or minhashed k-mers" << endl
			<< "    " << MODE_ALL_2_ALL << " - counting common k-mers - all samples in the database" << endl
			<< "    " << MODE_ALL_2_ALL_SPARSE << " - counting common k-mers - all samples in the database (sparse computation)" << endl
			<< "    " << MODE_ALL_2_ALL_PARTS << " - counting common k-mers - all samples in the database parts (sparse computation)" << endl
			//	<< "    " << MODE_DB_2_DB_SP << " - counting common k-mers - all samples between the databases (sparse computation)" << endl
			<< "    " << MODE_NEW_2_ALL << " - counting common k-mers - set of new samples versus database" << endl
			<< "    " << MODE_ONE_2_ALL << " - counting common k-mers - single sample versus database" << endl
			<< "    " << MODE_DISTANCE << " - calculating similarities/distances" << endl
			<< "    " << MODE_MINHASH << " - storing minhashed k-mers" << endl
			<< "Common options:" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl << endl
			<< "The meaning of other options and positional arguments depends on the selected mode. For more information, run:" << endl
			<< "kmer-db <mode> -help" << endl << endl;
	}
}
