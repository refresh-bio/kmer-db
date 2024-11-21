#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include "conversion.h"

#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <mutex>
#include <sstream>

#define LOG_VERBOSE(msg) { if (Log::getInstance(Log::LEVEL_VERBOSE)) { Log::getInstance(Log::LEVEL_VERBOSE) << msg << std::flush; }}
#define LOG_DEBUG(msg) { if (Log::getInstance(Log::LEVEL_DEBUG)) { Log::getInstance(Log::LEVEL_DEBUG) << msg << std::flush; }}
#define LOG_NORMAL(msg) { Log::getInstance(Log::LEVEL_NORMAL) << msg << std::flush; }


class LockedStream {
	std::ostream* out{ nullptr };
	std::recursive_mutex* mtx{ nullptr };

public:
//	LockedStream() : mtx() {}
	LockedStream() = default;
	LockedStream(std::ostream& out, std::recursive_mutex& mtx) : out(&out), mtx(&mtx) {}
	~LockedStream() {
		if (out) {
			out->flush();
			mtx->unlock();
		}
	}

	template <class T>
//	LockedStream& operator<< (const T& v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (bool v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (long v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (unsigned long v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (long long v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (unsigned long long v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (float v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (double v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (long double v)							{ if (out) { *out << v; }; return *this; }
//	LockedStream& operator<< (const void *v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (short v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (unsigned short v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (int v)							{ if (out) { *out << v; }; return *this; }
	LockedStream& operator<< (unsigned int v)							{ if (out) { *out << v; }; return *this; }
	
	LockedStream& operator<< (std::string v)							{ if (out) { *out << v; }; return *this; }
	
	LockedStream& operator<< (std::ostream& (*pf)(std::ostream&))	{ if (out) { *out << pf; }; return *this; }
	LockedStream& operator<< (std::ios& (*pf)(std::ios&))			{ if (out) { *out << pf; }; return *this; }
	LockedStream& operator<< (std::ios& (*pf)(std::ios_base&))		{ if (out) { *out << pf; }; return *this; }

	friend LockedStream& operator<<(LockedStream& os, const char* s);
	friend LockedStream& operator<<(LockedStream& os, const signed char* s);
	friend LockedStream& operator<<(LockedStream& os, const unsigned char* s);
};

inline LockedStream& operator<<(LockedStream& os, const char*s) { if (os.out) { *(os.out) << s; }; return os; }
inline LockedStream& operator<<(LockedStream& os, const signed char*s) { if (os.out) { *(os.out) << s; }; return os; }
inline LockedStream& operator<<(LockedStream& os, const unsigned char*s) { if (os.out) { *(os.out) << s; }; return os; }

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
	static Log& getInstance(int level) {
		static std::vector<std::shared_ptr<Log>> logs{
			std::shared_ptr<Log>(new Log()),
			std::shared_ptr<Log>(new Log()),
			std::shared_ptr<Log>(new Log())
		};
		
		return *logs[level];
	}

	// *****************************************************************************************
	template <class T>
	LockedStream operator<<(const T& v) {
		if (enabled) { 
			mtx.lock();
			out << v;
			return LockedStream(out, mtx);
		}
		return LockedStream();
	}

	// *****************************************************************************************
	LockedStream operator<< (std::ostream& (*pf)(std::ostream&)) {
		if (enabled) {
			mtx.lock();
			out << pf;
			return LockedStream(out, mtx);
		}
		return LockedStream();
	}

	// *****************************************************************************************
	LockedStream operator<< (std::ios& (*pf)(std::ios&)) {
		if (enabled) {
			mtx.lock();
			out << pf;
			return LockedStream(out, mtx);
		}
		return LockedStream();
	}

	// *****************************************************************************************
	LockedStream operator<< (std::ios& (*pf)(std::ios_base&)) {
		if (enabled) {
			mtx.lock();
			out << pf;
			return LockedStream(out, mtx);
		}
		return LockedStream();
	}

	static std::string formatLargeNumber(uint64_t num, int minWidth = 0);

	operator bool() const
	{
		return enabled;
	}

protected:
	bool enabled;
	std::ostream& out{std::cerr};
	std::recursive_mutex mtx;

	Log() : enabled(false) {}
};

