#pragma once

class Console
{
public:
	int parse(int argc, char** argv);




protected:
	int runBuildDatabase(const std::string& multipleKmcSamples, const std::string dbFilename);
	int runAllVsAll(const std::string& dbFilename, const std::string& similarityFile);
	int runOneVsAll(const std::string& dbFilename, const std::string& singleKmcSample, const std::string& similarityFile);
	int runListPatterns(const std::string& dbFilename, const std::string& patternFile);

	void showInstructions();

};