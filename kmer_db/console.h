#pragma once

class Console
{
public:
	int parse(int argc, char** argv);




protected:
	int numThreads;

	int runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile);
	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);
	int runMinHash(const std::string& multipleKmcSamples, float fraction);

	void showInstructions();

};