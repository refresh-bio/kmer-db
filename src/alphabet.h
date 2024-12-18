#pragma once
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <sstream>

enum class AlphabetType {
	nt,
	nt_preserve,
	aa,
	aa11_diamond,
	aa12_mmseqs,
	aa6_dayhoff,
	unknown
};

// *****************************************************************************************
//
class Alphabet {

public:
	struct Description {
		AlphabetType type;
		std::string name;
		std::string groups;
		bool preserveStrand;
	};

	Alphabet(const Description& desc) :
		type(desc.type),
		name(desc.name),
		preserveStrand(desc.preserveStrand),
		size(std::count(desc.groups.begin(), desc.groups.end(), ',') + 1),
		bitsPerSymbol(std::round(std::ceil(std::log2(size)))),
		maxKmerLen(64 / bitsPerSymbol - 1)  // highest bit in a value is reserved for the hash table use
		
	{ 
		std::fill_n(mapping, 256, -1);
		
		// translate comma-separated string into vector of strings	
		std::string line;
		std::vector<std::string> vec;
		std::stringstream ss(desc.groups);
		while (std::getline(ss, line, ',')) {
			vec.push_back(line);
		}
	
		for (int gi = 0; gi < vec.size(); ++gi) {
			const auto& group = vec[gi];

			for (unsigned char c : group) {
				mapping[std::tolower(c)] = mapping[std::toupper(c)] = gi;
			}
		}
	}

public:
	const AlphabetType type;
	const std::string name;
	const bool preserveStrand;
	const int size;
	const int bitsPerSymbol;
	const int maxKmerLen;
	
	int8_t map(char c) const { return mapping[c]; }

protected:
	int8_t mapping[256];
};

// *****************************************************************************************
//
class AlphabetFactory {
private:
	
	std::vector<Alphabet::Description> descriptions {
		{ AlphabetType::nt,				"nt",			"A,C,G,TU",									false },
		{ AlphabetType::nt_preserve,	"nt-preserve",	"A,C,G,TU",									true },
		{ AlphabetType::aa,				"aa",			"K,R,E,D,Q,N,C,G,H,I,L,V,M,F,Y,W,P,S,T,A",	true },
		{ AlphabetType::aa11_diamond,	"aa11_diamond",	"KREDQN,C,G,H,ILV,M,F,Y,W,P,STA",			true },
		{ AlphabetType::aa12_mmseqs,	"aa12_mmseqs",	"AST,C,DN,EQ,FY,G,H,IV,KR,LM,P,W",			true },
		{ AlphabetType::aa6_dayhoff,	"aa6_dayhoff",	"STPAG,NDEQ,HRK,MILV,FYW,C",					true },
	};
		
	AlphabetFactory() {}

public:

	// Creates an alphabet from enumeration
	Alphabet* create(AlphabetType type) {

		auto it = std::find_if(descriptions.begin(), descriptions.end(),
			[type](const Alphabet::Description& d) {return d.type == type; });

		if (it == descriptions.end()) {
			throw std::runtime_error("Invalid alphabet type");
		}

		return new Alphabet(*it);
	}

	// Creates an alphabet from string
	Alphabet* create(const std::string& name) {
		return create(str2type(name));
	}

	// Converts string to enumeration
	AlphabetType str2type(const std::string& name) {
		
		auto it = std::find_if(descriptions.begin(), descriptions.end(),
			[&name](const Alphabet::Description& d) {return d.name == name; });

		if (it == descriptions.end()) {
			throw std::runtime_error("Invalid alphabet type");
		}
		
		return it->type;
	}

	static AlphabetFactory& instance() {
		static AlphabetFactory factory;
		return factory;
	}
};