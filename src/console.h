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
template <class T>
class Param {
	bool set;
	T value;

public:
	Param(T def) : value(def) {}
	void setValue(T value) { this->value = value; set = true; }
	bool isSet() { return set; }
	T getValue() { return value; }
};


class Console
{
public:
	int parse(int argc, char** argv);

protected:
	int numThreads;
	int cacheBufferMb;

	int runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename, 
		InputFile::Format inputFormat, double filterValue, uint32_t kmerLength);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFilename);
	int runMinHash(const std::string& multipleKmcSamples, double fraction);
	int runDistanceCalculation(const std::string& similarityFilename);

	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);
	int runAnalyzeDatabase(const std::string& multipleKmcSamples, const std::string& dbFilename);

	void showInstructions();
};