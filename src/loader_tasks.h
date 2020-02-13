#pragma once


class InputFile;

// *****************************************************************************************
//
struct InputTask {
	size_t fileId;
	const std::string& filePath;
	std::shared_ptr<InputFile> file;

	
	InputTask(size_t fileId, const std::string& filePath) :
		fileId(fileId), filePath(filePath), file(nullptr) {
	}
};

// *****************************************************************************************
//
struct SampleTask {
	size_t fileId;
	const std::string& filePath;
	std::string sampleName;
	kmer_t *kmers;
	size_t kmersCount;
	uint32_t kmerLength;
	double fraction;
	int bufferId;

	SampleTask(size_t fileId, const std::string& filePath, const std::string& sampleName, int bufferId) :
		fileId(fileId), filePath(filePath), sampleName(sampleName), bufferId(bufferId) {}

};
