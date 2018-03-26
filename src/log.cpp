#include "log.h"
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include <sstream>
#include <iomanip>

using namespace std;

const int Log::LEVEL_DEBUG = 0;
const int Log::LEVEL_VERBOSE = 1;
const int Log::LEVEL_NORMAL = 2;

// *****************************************************************************************
//
Log::Log()
{
	enabled = false;
	out = &std::cerr;
}

// *****************************************************************************************
//
Log& Log::operator<< (std::ostream& (*pf)(std::ostream&))
{
	if (enabled) { 
		*this->out << pf; 
		out->flush();
	}
	return *this;
}

// *****************************************************************************************
//
Log& Log::operator<< (std::ios& (*pf)(std::ios&))
{
	if (enabled) { 
		*this->out << pf; 
		out->flush();
	}
	
	return *this;
}

// *****************************************************************************************
//
Log& Log::operator<< (std::ios_base& (*pf)(std::ios_base&))
{
	if (enabled) { 
		*this->out << pf;  
		out->flush();
	}

	return *this;
}

// *****************************************************************************************
//
std::string Log::formatLargeNumber(uint64_t num) {
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

	return out;
}