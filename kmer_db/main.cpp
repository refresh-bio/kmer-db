/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

Version: 1.0
Date   : 2018-02-10
*/

#include "kmer_db.h"
#include "log.h"

#include "console.h"

// *****************************************************************************************
//
int main(int argc, char **argv)
{
	Console console;
	console.parse(argc, argv);

	return 0;
}
