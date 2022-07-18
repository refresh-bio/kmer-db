/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include "kmer_db.h"
#include "log.h"

#include "console.h"

#include <stdexcept>

int main(int argc, char **argv)
{
	try {
		Console console;
		console.parse(argc, argv);
	}
	catch (std::runtime_error& err) {
		cout << "ERROR: " << err.what() << endl;
	}

	return 0;
}
