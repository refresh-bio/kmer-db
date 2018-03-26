#pragma once
#include "kmc_file_wrapper.h"
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

// *****************************************************************************************
//
class Console
{
public:
	int parse(int argc, char** argv);

protected:
	int numThreads;
	int cacheBufferMb;

	int runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename, KmcFileWrapper::Format inputFormat, double filterValue);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFilename);
	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);
	int runMinHash(const std::string& multipleKmcSamples, double fraction);
	int runDistanceCalculation(const std::string& similarityFilename);

	void showInstructions();
};