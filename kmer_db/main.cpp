#include "kmer_db.h"
#include "log.h"

#include "console.h"


/*
Example command line

minhash E:\Data\kmc250.list 0.1
build E:\Data\kmc250.list d:\kmer.db
build-minhash E:\Data\kmc250.list d:\kmer.db
all2all d:\kmer.db d:\matrix.csv
one2all d:\kmer.db E:\Data\kmc250\GCF_000171975.1_ASM17197v1_genomic_s d:\vector.csv
list-patterns d:\kmer.db d:\patterns.txt
distance d:\matrix.csv

*/


int main(int argc, char **argv)
{
	Console console;
	console.parse(argc, argv);
	

	return 0;
}

