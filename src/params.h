#pragma once

#include "input_file.h"
#include "sparse_filters.h"

#include <string>
#include <algorithm>
#include <vector>

// *****************************************************************************************
//

class Params {
public:


	
	enum Mode {
		build,
		minhash,
		all2all,
		all2all_sparse,
		all2all_parts,
		new2all,
		one2all,
		distance,
		unknown
	};
	
	const std::string MODE_BUILD = "build";
	const std::string MODE_MINHASH = "minhash";
	const std::string MODE_ALL_2_ALL = "all2all";
	const std::string MODE_ALL_2_ALL_SPARSE = "all2all-sp";
	const std::string MODE_ALL_2_ALL_PARTS = "all2all-parts";
	const std::string MODE_NEW_2_ALL = "new2all";
	const std::string MODE_ONE_2_ALL = "one2all";
	const std::string MODE_DISTANCE = "distance";

	const std::string SWITCH_HELP = "-help";
	const std::string SWITCH_KMC_SAMPLES = "-from-kmers";
	const std::string SWITCH_MINHASH_SAMPLES = "-from-minhash";
	const std::string SWITCH_MULTISAMPLE_FASTA = "-multisample-fasta";
	const std::string SWITCH_PHYLIP_OUT = "-phylip-out";
	const std::string SWITCH_EXTEND_DB = "-extend";
	const std::string SWITCH_SPARSE = "-sparse";
	const std::string SWITCH_NON_CANONICAL = "-non-canonical";

	const std::string OPTION_FRACTION = "-f";
	const std::string OPTION_FRACTION_START = "-f-start";
	const std::string OPTION_LENGTH = "-k";
	const std::string OPTION_VERBOSE = "-v";
	const std::string OPTION_DEBUG = "-vv";
	const std::string OPTION_THREADS = "-t";
	const std::string OPTION_READER_THREADS = "-rt";
	const std::string OPTION_BUFFER = "-buffer";

	const std::string OPTION_MAX = "-max";
	const std::string OPTION_MIN = "-min";

	const std::string OPTION_SAMPLE_ROWS = "-sample-rows";

	std::map<std::string, metric_fun_t> availableMetrics;

private:


	
public:
	double fraction{ 1.0 };
	uint32_t kmerLength{ 18 };

	int numThreads{ 0 };
	int numReaderThreads{ 0 };
	int cacheBufferMb{ 8 };
	bool multisampleFasta{ false };
	
	double fractionStart{ 0.0 };
	bool fractionSpecified{ false };
	
	bool sparseOut{ false };
	bool extendDb{ false };
	bool phylipOut{ false };
	bool nonCanonical{ false };

	int samplingSize{ 0 };
	metric_fun_t samplingCriterion;


	InputFile::Format inputFormat { InputFile::GENOME };
	Mode mode;

	std::vector<std::string> files;

	std::map<std::string, MetricFilter> metricFilters;
	KmerFilter kmerFilter;

	std::string metricName;

	
	Params();

	bool parse(int argc, char** argv);

	
	void showInstructions(Mode mode) const;

	
	bool findSwitch(std::vector<std::string>& params, const std::string& name) const {
		auto it = std::find(params.begin(), params.end(), name); // verbose mode
		if (it != params.end()) {
			params.erase(it);
			return true;
		}

		return false;
	}

	template <typename T>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& v) const {
		auto prevToEnd = std::prev(params.end());
		auto it = std::find(params.begin(), prevToEnd, name); // verbose mode
		if (it != prevToEnd) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				params.erase(it, it + 2);
				return true;
			}
		}

		return false;
	}

	template <typename T, typename U>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& value1, U& value2) {
		if (params.size() < 3) {
			return false;
		}

		auto stop = std::prev(params.end(), 2);
		auto it = find(params.begin(), stop, name); // verbose mode
		if (it != stop) {
			if (std::istringstream(*std::next(it)) >> value1
				&& std::istringstream(*std::next(it, 2)) >> value2) {
				params.erase(it, it + 3);
				return true;
			}
		}
		return false;
	}

	Mode str2mode(const std::string& str);

private:
	void parseFilters(std::vector<std::string>& params);

	bool parse_build(std::vector<std::string>& params);
	bool parse_all2all(std::vector<std::string>& params);
	bool parse_new2all(std::vector<std::string>& params);
	bool parse_distance(std::vector<std::string>& params);
	bool parse_minhash(std::vector<std::string>& params);
};
