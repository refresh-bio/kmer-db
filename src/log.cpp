#include "log.h"
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <sstream>
#include <iomanip>

using namespace std;

const int Log::LEVEL_DEBUG = 0;
const int Log::LEVEL_VERBOSE = 1;
const int Log::LEVEL_NORMAL = 2;


// ************************************************************************************
// NumericConversions statics
char NumericConversions::digits[];
NumericConversions::_si NumericConversions::_init;
uint64_t NumericConversions::powers10[];



// *****************************************************************************************
//
std::string Log::formatLargeNumber(uint64_t num, int minWidth) {
	std::string ret = "";

	do {
		uint64_t part = num % 1000uLL;
		num = num / 1000uLL;

		if (num > 0) {
			std::ostringstream oss;
			oss << "," << std::setw(3) << std::setfill('0') << part;
			ret = oss.str() + ret;
/*			auto s = std::to_string(part);
			if (s.length() < 3)
				ret = "," + std::string(3 - s.length(), '0') + s + ret;
			else
				ret = "," + s + ret;*/
		}
		else {
			ret = std::to_string(part) + ret;
		}

	} while (num > 0);

	int initialSpaces = minWidth - (int)ret.length();

	if (initialSpaces > 0) {
		ret = string(initialSpaces, ' ') + ret;
	}

	return ret;
}
