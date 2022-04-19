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
	std::string out = "";

	do {
		uint64_t part = num % 1000LL;
		num = num / 1000LL;

		if (num > 0) {
			std::ostringstream oss;
			oss << "," << std::setw(3) << std::setfill('0') << part;
			out = oss.str() + out;
		}
		else {
			out = std::to_string(part) + out;
		}

	} while (num > 0);

	int initialSpaces = minWidth - out.length();

	if (initialSpaces > 0) {
		out = string(initialSpaces, ' ') + out;
	}

	return out;
}
