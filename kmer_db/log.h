#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include <iostream>
#include <vector>
#include <memory>

#define LOG_NORMAL Log::getInstance(Log::LEVEL_NORMAL)
#define LOG_VERBOSE Log::getInstance(Log::LEVEL_VERBOSE)
#define LOG_DEBUG Log::getInstance(Log::LEVEL_DEBUG)

// *****************************************************************************************
//
class Log
{
public:
	static const int LEVEL_NORMAL;
	static const int LEVEL_VERBOSE;
	static const int LEVEL_DEBUG;

	void enable()	{ enabled = true; }
	void disable()	{ enabled = false; }
	
	// *****************************************************************************************
	//
	static Log& getInstance(int level) {
		static std::vector<std::shared_ptr<Log>> logs;
		if (logs.size() == 0) {
			logs.push_back(std::shared_ptr<Log>(new Log()));
			logs.push_back(std::shared_ptr<Log>(new Log()));
			logs.push_back(std::shared_ptr<Log>(new Log()));
		}

		return *logs[level];
	}

	// *****************************************************************************************
	//
	template <class T>
	Log& operator<<(T v) {
		if (enabled) { *out << v; }
		return *this;
	}

	Log& operator<< (std::ostream& (*pf)(std::ostream&));
	Log& operator<< (std::ios& (*pf)(std::ios&));
	Log& operator<< (std::ios_base& (*pf)(std::ios_base&));

	static std::string formatLargeNumber(uint64_t num);

protected:
	bool enabled;
	std::ostream* out;

	Log();
};



