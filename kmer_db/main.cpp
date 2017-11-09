// kmer_db.cpp : Defines the entry point for the console application.
//



using namespace std;


#include "kmer_db.h"
#include "tests.h"
#include "log.h"

#include "console.h"


int main(int argc, char **argv)
{
	//Log::getInstance(Log::LEVEL_DEBUG).enable();
	Console console;
	console.parse(argc, argv);
	

	return 0;
}

