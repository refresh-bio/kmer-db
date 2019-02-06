/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#pragma once
#include "kmc_file_wrapper.h"
#include <stdexcept>

// *****************************************************************************************
//

class Params {
public:


	static const string MODE_BUILD;
	static const string MODE_MINHASH;
	static const string MODE_ALL_2_ALL;
	static const string MODE_NEW_2_ALL;
	static const string MODE_ONE_2_ALL;
	static const string MODE_DISTANCE;

	static const string SWITCH_KMC_SAMPLES;
	static const string SWITCH_MINHASH_SAMPLES;
	
	static const string OPTION_FILTER;
	static const string OPTION_LENGTH;
	static const string OPTION_VERBOSE;
	static const string OPTION_THREADS;
	static const string OPTION_BUFFER;
};

class Console
{
public:
	int parse(int argc, char** argv);

protected:
	int numThreads;
	int cacheBufferMb;

	int runBuildDatabase(const std::string& multipleSamples, const std::string dbFilename, 
		InputFile::Format inputFormat, double filterValue, uint32_t kmerLength);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runNewVsAll(const std::string& dbFilename, const std::string& multipleSamples, const std::string& similarityFilename, InputFile::Format inputFormat);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFilename, InputFile::Format inputFormat);

	int runMinHash(const std::string& multipleKmcSamples, double fraction);
	int runDistanceCalculation(const std::string& similarityFilename);

	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);
	int runAnalyzeDatabase(const std::string& multipleKmcSamples, const std::string& dbFilename);

	void showInstructions();

	bool findSwitch(std::vector<std::string>& params, const std::string& name) {
		auto it = find(params.begin(), params.end(), name); // verbose mode
		if (it != params.end()) {
			params.erase(it);
			return true;
		}

		return false;
	}

	template <typename T>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& v) {
		auto prevToEnd = std::prev(params.end());
		auto it = find(params.begin(), prevToEnd, name); // verbose mode
		if (it != prevToEnd) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				params.erase(it, it + 2);
				return true;
			}
		}

		return false;
	}
};