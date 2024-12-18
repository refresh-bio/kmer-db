#include "params.h"
#include "types.h"
#include "version.h"

#include <iostream>
#include <stdexcept>
#include <cstdint>

using namespace std;


// *****************************************************************************************
//
Params::Params() {
	availableMetrics["jaccard"] = [](num_kmers_t common, num_kmers_t cnt1, num_kmers_t cnt2, int kmerLength) -> double { return (double)common / (cnt1 + cnt2 - common); };
	availableMetrics["min"] = [](num_kmers_t common, num_kmers_t cnt1, num_kmers_t cnt2, int kmerLength) -> double { return (double)common / std::min(cnt1, cnt2); };
	availableMetrics["max"] = [](num_kmers_t common, num_kmers_t cnt1, num_kmers_t cnt2, int kmerLength) -> double { return (double)common / std::max(cnt1, cnt2); };
	availableMetrics["cosine"] = [](num_kmers_t common, num_kmers_t cnt1, num_kmers_t cnt2, int kmerLength) -> double { return (double)common / sqrt(cnt1 * cnt2); };
	availableMetrics["mash"] = [](num_kmers_t common, num_kmers_t queryCnt, num_kmers_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / (queryCnt + dbCnt - common);
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
	};
	availableMetrics["ani"] = [](num_kmers_t common, num_kmers_t queryCnt, num_kmers_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / (queryCnt + dbCnt - common);
		double d_mash = (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
		return 1.0 - d_mash;
	};
	availableMetrics["ani-shorter"] = [](num_kmers_t common, num_kmers_t queryCnt, num_kmers_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / min(queryCnt, dbCnt);
		double d_mash = (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
		return 1.0 - d_mash;
	};

	availableMetrics["mash-query"] = [](num_kmers_t common, num_kmers_t queryCnt, num_kmers_t dbCnt, int kmerLength) -> double {
		double d_jaccard = (double)common / queryCnt;
		return  (d_jaccard == 0) ? 1.0 : (-1.0 / kmerLength) * log((2 * d_jaccard) / (d_jaccard + 1));
	};

	availableMetrics["num-kmers"] = [](num_kmers_t common, num_kmers_t queryCnt, num_kmers_t dbCnt, int kmerLength) -> double {
		return (double)common;
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

	if (findSwitch(params, SWITCH_VERSION)) {
		LOG_NORMAL(VERSION);
		return false;
	}
	
	showHeader();

	bool helpWanted = findSwitch(params, SWITCH_HELP);
	if (params.size() == 0) {
		// no arguments or -help switch already consumed
		showInstructions(Mode::unknown);
		return false;
	}

	mode = str2mode(params.front());
	params.erase(params.begin());
	
	if (helpWanted || params.size() == 0 || mode == Mode::unknown) {
		// help for particular mode
		showInstructions(mode);
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

		// parse mode-specific params
		switch (mode) {
		case Mode::build :
			parse_build(params); break;
		case Mode::all2all:
		case Mode::all2all_sparse:
		case Mode::all2all_parts:
			parse_all2all(params); break;
		case Mode::new2all:
		case Mode::one2all:
			parse_new2all(params); break;
		case Mode::distance:
			parse_distance(params); break;
		case Mode::minhash:
			parse_minhash(params); break;
		default:
			showInstructions(mode);
			return false;
		}

		// set default fration in minhash mode if not specified
		if (this->mode == Mode::minhash && !fractionSpecified) {
			fraction = 0.01;
		}

		files.insert(files.begin(), params.begin(), params.end());
	}

	return true;
}

// *****************************************************************************************
//
void Params::showHeader() const {

	LOG_NORMAL("Kmer-db version " << VERSION << " (" << DATE << ")" << endl
		<< "S. Deorowicz, A. Gudys, M. Dlugosz, M. Kokot, and A. Danek (c) 2018" << endl << endl);
}

// *****************************************************************************************
//
void Params::showInstructions(Mode mode) const {

	if (mode == Mode::build) {
		LOG_NORMAL(
			   "Building a database:" << endl
	 << "    kmer-db " << MODE_BUILD
			<< " [" << OPTION_LENGTH << " <kmer-length>]"
			<< " [" << OPTION_FRACTION << " <fraction>]"
			<< " [" << SWITCH_MULTISAMPLE_FASTA << "]"
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_ALPHABET << "<alphabet_type>]"
			<< " [" << SWITCH_PRESERVE_STRAND << "]"
			<< " [" << OPTION_THREADS << " <threads>] <samples> <database>" << endl

			<< "    kmer-db " << MODE_BUILD << " " << SWITCH_KMC_SAMPLES
			<< " [" << OPTION_FRACTION << " <fraction>]"
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_THREADS << " <threads>] <samples> <database>" << endl

			<< "    kmer-db " << MODE_BUILD << " " << SWITCH_MINHASH_SAMPLES
			<< " [" << SWITCH_EXTEND_DB << "]"
			<< " [" << OPTION_THREADS << " <threads>] <samples> <database>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    samples (input) - one of the following:" << endl 
			<< "        a. FASTA file (fa, fna, fasta, fa.gz, fna.gz, fasta.gz) with one or multiple (" << SWITCH_MULTISAMPLE_FASTA << ") samples" << endl
			<< "        b. file with list of samples in one of the formats: " << endl
			<< "            * FASTA genomes/reads (default)," << endl 
			<< "            * KMC k-mers (" << SWITCH_KMC_SAMPLES << ")," << endl
			<< "            * minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << ")" << endl
			<< "    database (output) - file with generated k-mer database," << endl

			<< "Options: " << endl
			<< "    " << OPTION_LENGTH << " <kmer_length> - length of k-mers (default: 18, maximum depends on the alphabet - 31 for default nt)" << endl
			<< "    " << OPTION_FRACTION << " <fraction> - fraction of all k-mers to be accepted by the minhash filter (default: 1)" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl
			<< "    " << SWITCH_EXTEND_DB << " - extend the existing database with new samples" << endl
			<< "    " << OPTION_ALPHABET << " - alphabet:" << endl
			<< "        * nt (4 symbol nucleotide with indistinguishable T/U; default)" << endl 
			<< "        * aa (20 symbol amino acid)" << endl
			<< "        * aa12_mmseqs (amino acid reduced to 12 symbols as in MMseqs: AST,C,DN,EQ,FY,G,H,IV,KR,LM,P,W" << endl
			<< "        * aa11_diamond (amino acid reduced to 11 symbols as in Diamond: KREDQN,C,G,H,ILV,M,F,Y,W,P,STA" << endl
			<< "        * aa6_dayhoff (amino acid reduced to 6 symbols as proposed by Dayhoff: STPAG,NDEQ,HRK,MILV,FYW,C" << endl
			<< "    " << SWITCH_PRESERVE_STRAND << " - preserve strand instead of taking canonical k-mers (allowed only in nt alphabet; default: off)" << endl 
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl << endl);
	}
	else if (mode == Mode::all2all) {
		LOG_NORMAL(
			   "Counting common k-mers for all the samples in the database:" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << SWITCH_SPARSE << " [" << OPTION_MIN << " [<criterion>:]<value>]* [" << OPTION_MAX << " [<criterion>:]<value>]* ]"
			<< " <database> <common_table>" << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << SWITCH_SPARSE << " - produce sparse matrix as output" << endl
			<< "    " << OPTION_MIN << " [<criterion>:]<value> - retains elements with <criterion> greater than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_MAX << " [<criterion>:]<value> - retains elements with <criterion> lower than or equal to <value> (see details below)" << endl
			<< "<criterion> can be num-kmers (number of common k-mers) or one of the distance/similarity measures: jaccard, min, max, cosine, mash, ani, ani-shorder." << endl
			<< "No <criterion> indicates num-kmers. Multiple filters can be combined." << endl
			<< endl);
	}
	else if (mode == Mode::all2all_sparse) {
		LOG_NORMAL(
			   "Counting common k-mers for all the samples in the database (sparse computation):" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL_SPARSE
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << OPTION_MIN << " [<criterion>:]<value>]* [" << OPTION_MAX << " [<criterion>:]<value>]* " << endl
			<< " <database> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    common_table (output) - comma-separated table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << OPTION_MIN << " [<criterion>:]<value> - retains elements with <criterion> greater than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_MAX << " [<criterion>:]<value> - retains elements with <criterion> lower than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_SAMPLE_ROWS << " [<criterion>:]<count> - retains <count> elements in every row using one of the strategies:" << endl 
			<< "        (i) random selection (no <criterion>)," << endl 
			<< "       (ii) the best elements with respect to <criterion>" << endl
			<< "<criterion> can be num-kmers (number of common k-mers) or one of the distance/similarity measures: jaccard, min, max, cosine, mash, ani, ani-shorder." << endl
			<< "No <criterion> indicates num-kmers (filtering) or random elements selection (sampling). Multiple filters can be combined." << endl
			<< endl);
	}
	else if (mode == Mode::all2all_parts) {
		LOG_NORMAL(
			   "Counting common k-mers for all the samples in the database parts (sparse computation):" << endl
			<< "    kmer-db " << MODE_ALL_2_ALL_PARTS
			<< " [" << OPTION_BUFFER << " <size_mb>]"
			<< " [" << OPTION_THREADS << " <threads>]" 
			<< " [" << OPTION_MIN << " [<criterion>:]<value>]* [" << OPTION_MAX << " [<criterion>:]<value>]* "
			<< " <db_list> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    db_list (input) - file containing list of database file names" << endl
			<< "    common_table (output) - CSV table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << OPTION_BUFFER << " <size_mb> - size of cache buffer in megabytes" << endl
			<< "                      (use L3 size for Intel CPUs and L2 for AMD to maximize performance; default: 8)" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << OPTION_MIN << " [<criterion>:]<value> - retains elements with <criterion> greater than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_MAX << " [<criterion>:]<value> - retains elements with <criterion> lower than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_SAMPLE_ROWS << " [<criterion>:]<count> - retains <count> elements in every row using one of the strategies:" << endl 
			<< "        (i) random selection (no <criterion>)," << endl
			<< "       (ii) the best elements with respect to <criterion>" << endl
			<< "<criterion> can be num-kmers (number of common k-mers) or one of the distance/similarity measures: jaccard, min, max, cosine, mash, ani, ani-shorder." << endl
			<< "No <criterion> indicates num-kmers (filtering) or random elements selection (sampling). Multiple filters can be combined." << endl
			<< endl);
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
		LOG_NORMAL(
			   "Counting common kmers between set of new samples and all the samples in the database:" << endl
			<< "    kmer-db " << MODE_NEW_2_ALL
			<< " [" << SWITCH_MULTISAMPLE_FASTA << " | " << SWITCH_KMC_SAMPLES << " | " << SWITCH_MINHASH_SAMPLES << "]"
			<< " [" << SWITCH_SPARSE << " [" << OPTION_MIN << " [<criterion>:]<value>]* [" << OPTION_MAX << " [<criterion>:]<value>]* ] " 
			<< " [" << OPTION_THREADS << " <threads>] <database> <samples> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file" << endl
			<< "    samples (input) - one of the following : " << endl 
			<< "        a. FASTA file (fa, fna, fasta, fa.gz, fna.gz, fasta.gz) with one or multiple (" << SWITCH_MULTISAMPLE_FASTA << ") samples" << endl
			<< "        b. file with list of samples in one of the formats: " << endl
			<< "            * FASTA genomes/reads (default)," << endl
			<< "            * KMC k-mers (" << SWITCH_KMC_SAMPLES << ")," << endl
			<< "            * minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << ")" << endl
			<< "    common_table (output) - CSV table with number of common k-mers" << endl

			<< "Options:" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl
			<< "    " << SWITCH_SPARSE << " - outputs a sparse matrix" << endl
			<< "    " << OPTION_MIN << " [<criterion>:]<value> - retains elements with <criterion> greater than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_MAX << " [<criterion>:]<value> - retains elements with <criterion> lower than or equal to <value> (see details below)" << endl
			<< "<criterion> can be num-kmers (number of common k-mers) or one of the distance/similarity measures: jaccard, min, max, cosine, mash, ani, ani-shorder." << endl
			<< "No <criterion> indicates num-kmers. Multiple filters can be combined." << endl
			<< endl);
	}
	else if (mode == Mode::one2all) {
		LOG_NORMAL(
			   "Counting common kmers between single sample and all the samples in the database:" << endl
			<< "    kmer-db " << MODE_ONE_2_ALL
			<< " [" << SWITCH_KMC_SAMPLES << " | " << SWITCH_MINHASH_SAMPLES << "] <database> <sample> <common_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    database (input) - k-mer database file." << endl
			<< "    sample (input) - query sample in one of the supported formats:" << endl
			<< "                     FASTA genomes/reads (default), KMC k-mers (" << SWITCH_KMC_SAMPLES << "), or minhashed k-mers (" << SWITCH_MINHASH_SAMPLES << "), " << endl
			<< "    common_table (output) - CSV table with number of common k-mers." << endl << endl);
	}
	else if (mode == Mode::distance) {
		LOG_NORMAL(
			   "Calculating similarities/distances on the basis of common k-mers:" << endl
			<< "    kmer-db " << MODE_DISTANCE << " <measure>"
			<< " [" << SWITCH_SPARSE << " [" << OPTION_MIN << "[<criterion>:]<value>]* [" << OPTION_MAX << " [<criterion>:]<value>]* ] "
			<< "<common_table> <output_table>" << endl << endl

			<< "Positional arguments:" << endl
			<< "    measure - name of the similarity/distance measure to be calculated, one of the following:" << endl
			<< "               jaccard, min, max, cosine, mash, ani, ani-shorter." << endl
			<< "    common_table (input) - CSV table with a number of common k-mers" << endl
			<< "    output_table (output) - CSV table with calculated distances" << endl

			<< "Options:" << endl

			<< "    " << SWITCH_PHYLIP_OUT << " - store output distance matrix in a Phylip format" << endl
			<< "    " << SWITCH_SPARSE << " - outputs a sparse matrix (only for dense input matrices - sparse input always produce sparse output)" << endl
			<< "    " << OPTION_MIN << " [<criterion>:]<value> - retains elements with <criterion> greater than or equal to <value> (see details below)" << endl
			<< "    " << OPTION_MAX << " [<criterion>:]<value> - retains elements with <criterion> lower than or equal to <value> (see details below)" << endl
			<< "<criterion> can be num-kmers (number of common k-mers) or one of the distance/similarity measures: jaccard, min, max, cosine, mash, ani, ani-shorder." << endl
			<< "If no criterion is specified, measure argument is used by default. Multiple filters can be combined." << endl
			<< endl);

	}
	else if (mode == Mode::minhash) {
		LOG_NORMAL(
			   "Storing minhashed k-mers:" << endl
			<< "    kmer-db " << MODE_MINHASH 
			<< " [" << OPTION_FRACTION << " <fraction> ]" 
			<< " [" << OPTION_LENGTH << " <kmer_length>]"  
			<< " [" << SWITCH_MULTISAMPLE_FASTA << "]"
			<< " [" << OPTION_ALPHABET << "<alphabet_type>]"
			<< " [" << SWITCH_PRESERVE_STRAND << "]"
			<< " <samples>" << endl

			<< "    kmer-db " << MODE_MINHASH << " " << SWITCH_KMC_SAMPLES 
			<< " [" << OPTION_FRACTION << " <fraction>]"
			<< " <samples>" << endl << endl
			
			<< "Positional arguments:" << endl
			<< "    samples (input) - one of the following : " << endl 
			<< "        a. FASTA file (fa, fna, fasta, fa.gz, fna.gz, fasta.gz) with one or multiple (" << SWITCH_MULTISAMPLE_FASTA << ") samples" << endl
			<< "        b. file with list of samples in one of the formats: " << endl
			<< "            * FASTA genomes/reads (default)," << endl
			<< "            * KMC k-mers (" << SWITCH_KMC_SAMPLES << ")," << endl << endl

			<< "Options:" << endl
			<< "    " << OPTION_FRACTION << " <fraction> - fraction of all k-mers to be accepted by the minhash filter (default: 0.01)" << endl
			<< "    " << OPTION_LENGTH << " <kmer_length> - length of k-mers (default: 18, maximum: 30)" << endl
			<< "    " << SWITCH_MULTISAMPLE_FASTA << " - each sequence in a FASTA file is treated as a separate sample" << endl << endl
			<< "    " << OPTION_ALPHABET << " - alphabet:" << endl
			<< "        * nt (4 symbol nucleotide with indistinguishable T/U; default)" << endl
			<< "        * aa (20 symbol amino acid)" << endl
			<< "        * aa12_mmseqs (amino acid reduced to 12 symbols as in MMseqs: AST,C,DN,EQ,FY,G,H,IV,KR,LM,P,W" << endl
			<< "        * aa11_diamond (amino acid reduced to 11 symbols as in Diamond: KREDQN,C,G,H,ILV,M,F,Y,W,P,STA" << endl
			<< "        * aa6_dayhoff (amino acid reduced to 6 symbols as proposed by Dayhoff: STPAG,NDEQ,HRK,MILV,FYW,C" << endl
			<< "    " << SWITCH_PRESERVE_STRAND << " - preserve strand instead of taking canonical k-mers (allowed only in nt alphabet; default: off)" << endl << endl

			<< "For each sample from the list, a binary file with *.minhash* extension containing filtered k-mers is created" << endl << endl);
	}
	else {
		LOG_NORMAL(
			   "USAGE" << endl
			<< "    kmer-db <mode> [options] <positional arguments>" << endl << endl

			<< "Modes:" << endl
			<< "    " << MODE_BUILD << " - building a database from FASTA genomes/reads, k-mers, or minhashed k-mers" << endl
			<< "    " << MODE_ALL_2_ALL << " - counting common k-mers - all samples in the database" << endl
			<< "    " << MODE_ALL_2_ALL_SPARSE << " - counting common k-mers - all samples in the database (sparse computation)" << endl
			<< "    " << MODE_ALL_2_ALL_PARTS << " - counting common k-mers - all samples in the database parts (sparse computation)" << endl
			<< "    " << MODE_NEW_2_ALL << " - counting common k-mers - set of new samples versus database" << endl
			<< "    " << MODE_ONE_2_ALL << " - counting common k-mers - single sample versus database" << endl
			<< "    " << MODE_DISTANCE << " - calculating similarities/distances" << endl
			<< "    " << MODE_MINHASH << " - storing minhashed k-mers" << endl
			<< "Common options:" << endl
			<< "    " << OPTION_THREADS << " <threads> - number of threads (default: number of available cores)" << endl << endl
			<< "The meaning of other options and positional arguments depends on the selected mode. For more information, run:" << endl
			<< "kmer-db <mode> -help" << endl << endl);
	}
}


// *****************************************************************************************
//
void Params::parseFilters(std::vector<std::string>& params) {
	std::vector<string> options{ OPTION_MIN, OPTION_MAX };

	for (int i = 0; i < options.size(); ++i) {

		string value_str;
		while (findOption(params, options[i], value_str)) {
			string metric;
			double value;

			auto beg = value_str.begin();
			auto sep = value_str.rfind(':');
			if (sep != string::npos) {
				metric = string(value_str.begin(), value_str.begin() + sep);
				beg = value_str.begin() + sep + 1;
			}
			else {
				metric = "num-kmers";
			}

			istringstream iss(string(beg, value_str.end()));
			if (!(iss >> value)) {
				throw runtime_error("Filtering error - unable to parse numerical value: " + value_str);
			}

			if (metric == "num-kmers") {
				kmerFilter.bounds[i] = lrint(value);
			}
			else if (availableMetrics.contains(metric)) {
				metricFilters[metric].metric = availableMetrics[metric];
				metricFilters[metric].bounds[i] = value;
			}
			else {
				throw runtime_error("Filtering error - unknown metric: " + metric);
			}
		}
	}
}

// *****************************************************************************************
//
bool Params::parse_build(std::vector<string>& params) {

	bool switchKmc = findSwitch(params, SWITCH_KMC_SAMPLES);
	bool switchMinhash = findSwitch(params, SWITCH_MINHASH_SAMPLES);

	if (!switchMinhash) {
		fractionSpecified = findOption(params, OPTION_FRACTION, fraction);	// minhash fraction
		findOption(params, OPTION_FRACTION_START, fractionStart);			// minhash fraction start value

		if (!switchKmc) {
			multisampleFasta = findSwitch(params, SWITCH_MULTISAMPLE_FASTA);
			
			inputFormat = InputFile::GENOME;

			string alphabetName;
			if (findOption(params, OPTION_ALPHABET, alphabetName)) {
				alphabet = std::shared_ptr<Alphabet>(AlphabetFactory::instance().create(alphabetName));
			}

			// override alphabet type when needed
			if (findSwitch(params, SWITCH_PRESERVE_STRAND)) {
				if (alphabet->type == AlphabetType::nt) {
					alphabet = std::shared_ptr<Alphabet>(AlphabetFactory::instance().create(AlphabetType::nt_preserve));
				}
				else {
					throw std::runtime_error("Switch " + SWITCH_PRESERVE_STRAND + " applies only to nt alphabet");
				}
			}

			findOption(params, OPTION_LENGTH, kmerLength);				// kmer length
			if (kmerLength > alphabet->maxKmerLen) {
				throw std::runtime_error("K-mer length for the given alphabet cannot exceed " + alphabet->maxKmerLen);
			}
		}
		else {
			inputFormat = InputFile::KMC;
			kmerLength = 0;
		}
	}
	else {
		if (switchKmc) {
			throw std::runtime_error(
				SWITCH_KMC_SAMPLES + " and " + SWITCH_MINHASH_SAMPLES + " switches exclude one another.");
		}

		inputFormat = InputFile::MINHASH;
		fraction = 1.0;
		kmerLength = 0;
	}

	extendDb = findSwitch(params, SWITCH_EXTEND_DB);
	

	return true;
}

// *****************************************************************************************
//
bool Params::parse_all2all(std::vector<string>& params) {

	findOption(params, OPTION_BUFFER, cacheBufferMb);	// size of temporary buffer in megabytes
	if (cacheBufferMb <= 0) {
		cacheBufferMb = 8;
	}

	findOption(params, OPTION_BUBBLE_SIZE, bubbleSize);	// min size of a bubble

	sparseOut = findSwitch(params, SWITCH_SPARSE);

	// parse filters
	if (sparseOut || mode == Mode::all2all_parts || mode == all2all_sparse) {
		parseFilters(params);
	}

	if (mode == Mode::all2all_parts || mode == Mode::all2all_sparse) {
		string value_str;
		if (findOption(params, OPTION_SAMPLE_ROWS, value_str)) {

			auto beg = value_str.begin();
			auto sep = value_str.rfind(':');
			if (sep != string::npos) {
				string measure = string(value_str.begin(), value_str.begin() + sep);

				if (availableMetrics.contains(measure)) {
					samplingCriterion = availableMetrics.at(measure);
				}
				else {
					throw runtime_error("Sampling parameters error - unknown measure: " + measure);
				}

				beg = value_str.begin() + sep + 1;
			}

			istringstream iss(string(beg, value_str.end()));
			if (!(iss >> samplingSize)) {
				throw runtime_error("Sampling parameters error - unable to parse numerical value: " + value_str);
			}
		}
		

	}

	return true;
}

// *****************************************************************************************
//
bool Params::parse_new2all(std::vector<string>& params) {

	bool switchKmc = findSwitch(params, SWITCH_KMC_SAMPLES);
	bool switchMinhash = findSwitch(params, SWITCH_MINHASH_SAMPLES);

	if (!switchMinhash) {
		if (!switchKmc) {
			multisampleFasta = findSwitch(params, SWITCH_MULTISAMPLE_FASTA);
			inputFormat = InputFile::GENOME;
		}
		else {
			inputFormat = InputFile::KMC;
		}
	}
	else {
		if (switchKmc) {
			throw std::runtime_error(
				SWITCH_KMC_SAMPLES + " and " + SWITCH_MINHASH_SAMPLES + " switches exclude one another.");
		}

		inputFormat = InputFile::MINHASH;
	}

	if (mode == Mode::new2all) {
		sparseOut = findSwitch(params, SWITCH_SPARSE);

		// parse filters
		if (sparseOut) {
			parseFilters(params);
		}
	}

	return true;
}

// *****************************************************************************************
//
bool Params::parse_distance(std::vector<std::string>& params) {

	sparseOut = findSwitch(params, SWITCH_SPARSE);

	phylipOut = findSwitch(params, SWITCH_PHYLIP_OUT);
	if (phylipOut) {
		sparseOut = false;
	}

	std::vector<string> options{ OPTION_MIN, OPTION_MAX };

	for (int i = 0; i < options.size(); ++i) {

		string value_str;
		while (findOption(params, options[i], value_str)) {
			string metric;
			double value;

			auto beg = value_str.begin();
			auto sep = value_str.rfind(':');
			if (sep != string::npos) {
				metric = string(value_str.begin(), value_str.begin() + sep);
				beg = value_str.begin() + sep + 1;
			}
			else {
				metric = "?";
			}

			istringstream iss(string(beg, value_str.end()));
			if (!(iss >> value)) {
				throw runtime_error("Filtering error - unable to parse numerical value: " + value_str);
			}

			if (metric == "num-kmers") {
				kmerFilter.bounds[i] = lrint(value);
			}
			else if (availableMetrics.contains(metric)) {
				metricFilters[metric].metric = availableMetrics[metric];
				metricFilters[metric].bounds[i] = value;
			}
			else if (metric == "?") {
				// metric to be filled later
				metricFilters[metric].bounds[i] = value;
			}
			else {
				throw runtime_error("Filtering error - unknown metric: " + metric);
			}
		}
	}

	if (params.empty()) {
		throw runtime_error("No distance/similarity metric specified");
	}

	metricName = params.front();
	params.erase(params.begin());

	// replace ? with real metric
	auto it = metricFilters.find("?");
	if (it != metricFilters.end()) {
		MetricFilter mf(it->second);
		mf.metric = availableMetrics[metricName];
		metricFilters[metricName] = mf;
		metricFilters.erase(it);
	}

	return true;
}

// *****************************************************************************************
//
bool Params::parse_minhash(std::vector<string>& params) {

	fractionSpecified = findOption(params, OPTION_FRACTION, fraction);		// minhash fraction
	findOption(params, OPTION_FRACTION_START, fractionStart);				// minhash fraction start value

	if (findSwitch(params, SWITCH_KMC_SAMPLES)) {
		inputFormat = InputFile::KMC;
		kmerLength = 0;
	}
	else {
		multisampleFasta = findSwitch(params, SWITCH_MULTISAMPLE_FASTA);
		findOption(params, OPTION_LENGTH, kmerLength);				// kmer length

		inputFormat = InputFile::GENOME;

		string alphabetName;
		if (findOption(params, OPTION_ALPHABET, alphabetName)) {
			alphabet = std::shared_ptr<Alphabet>(AlphabetFactory::instance().create(alphabetName));
		}

		// override alphabet type when needed
		if (findSwitch(params, SWITCH_PRESERVE_STRAND)) {
			if (alphabet->type == AlphabetType::nt) {
				alphabet = std::shared_ptr<Alphabet>(AlphabetFactory::instance().create(AlphabetType::nt_preserve));
			}
			else {
				throw std::runtime_error("Switch " + SWITCH_PRESERVE_STRAND + " applies only to nt alphabet");
			}
		}

		if (kmerLength > alphabet->maxKmerLen) {
			throw std::runtime_error("K-mer length for the given alphabet cannot exceed " + alphabet->maxKmerLen);
		}
	}

	return true;
}