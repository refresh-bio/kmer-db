/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#pragma once
#include "params.h"
#include <stdexcept>
#include <functional>
#include <memory>

// *****************************************************************************************
class usage_error : public std::runtime_error {
	Params::Mode mode;
public:
	usage_error(Params::Mode mode) : std::runtime_error(""), mode(mode) {}
	Params::Mode getMode() const { return mode; }
};

// *****************************************************************************************
class Console
{
public:

	virtual void run(const Params& params) = 0;
	virtual ~Console() {}
};
	
// *****************************************************************************************
class BuildConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class All2AllConsole : public Console {
public:
	virtual void run(const Params& params) override;
};


class All2AllSparseConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class All2AllPartsConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class New2AllConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class One2AllConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class Db2DbConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class MinhashConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

class DistanceConsole : public Console {
public:
	virtual void run(const Params& params) override;
};

// *****************************************************************************************
class ConsoleFactory {
public:
	
	static std::unique_ptr<Console> create(Params::Mode mode) {

		Console* p = nullptr;

		if (mode == Params::Mode::build) {
			p = new BuildConsole();
		} else if (mode == Params::Mode::all2all) {
			p = new All2AllConsole();
		} else if (mode == Params::Mode::all2all_parts) {
			p = new All2AllPartsConsole();
		} else if (mode == Params::Mode::all2all_sparse) {
			p = new All2AllSparseConsole();
		} else if (mode == Params::Mode::new2all) {
			p = new New2AllConsole();
		} else if (mode == Params::Mode::one2all) {
			p = new One2AllConsole();
		} else if (mode == Params::Mode::distance) {
			p = new DistanceConsole();
		} else if (mode == Params::Mode::minhash) {
			p = new MinhashConsole();
		} 

		return std::unique_ptr<Console>(p);
	}
};