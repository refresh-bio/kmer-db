#pragma once
#include "genome_input_file.h"
#include "kmc_input_file.h"
#include "minhashed_input_file.h"

class InputFileFactory {
public:
	
	static InputFile* create(
		InputFile::Format format, 
		std::shared_ptr<AbstractFilter> filter,
		std::shared_ptr<Alphabet> alphabet) {
		
		if (format == InputFile::MINHASH) {
			return new MihashedInputFile();
		}
		else {
			std::shared_ptr<MinHashFilter> minhashFilter = std::dynamic_pointer_cast<MinHashFilter>(filter);
			std::shared_ptr<NullFilter> nullFilter = std::dynamic_pointer_cast<NullFilter>(filter);

			if (!minhashFilter) {
				throw std::runtime_error("Only MinHashFilter is currently supported");
			}
			else if (nullFilter) {
				
				if (format == InputFile::GENOME) {
					return new GenomeInputFile<NullFilter>(nullFilter, alphabet);
				}
				else if (format == InputFile::KMC) {
					return new KmcInputFile<NullFilter>(nullFilter);
				}
				else {
					throw std::runtime_error("Unsupported input type");
				}
			}
			else {
				if (format == InputFile::GENOME) {
					return new GenomeInputFile<MinHashFilter>(minhashFilter, alphabet);
				}
				else if (format == InputFile::KMC) {
					return new KmcInputFile<MinHashFilter>(minhashFilter);
				}
				else {
					throw std::runtime_error("Unsupported input type");
				}
			}
		}
	}
};