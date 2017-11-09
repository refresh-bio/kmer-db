#include "kmer_db.h"
#include "tests.h"
#include "log.h"

#include "console.h"


/*
Example command line


--build E:\Data\kmc250.list d:\kmer.db
--all2all d:\kmer.db d:\matrix.csv
--one2all d:\kmer.db E:\Data\kmc250\GCF_000171975.1_ASM17197v1_genomic_s d:\vector.csv

*/


int main(int argc, char **argv)
{
	//Log::getInstance(Log::LEVEL_DEBUG).enable();
	Console console;
	console.parse(argc, argv);
	

	return 0;
}

