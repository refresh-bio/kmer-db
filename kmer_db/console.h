#pragma once

class Console
{
public:
	int parse(int argc, char** argv);




protected:
	int numThreads;

	int runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename, bool loadMinhash);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFilename);
	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);
	int runMinHash(const std::string& multipleKmcSamples, float fraction);
	int runDistanceCalculation(const std::string& similarityFilename);

	void showInstructions();

};