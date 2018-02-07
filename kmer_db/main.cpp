#include "kmer_db.h"
#include "log.h"

#include "console.h"


/*
Example command lines

build E:\Data\kmc250.list d:\kmer.db
all2all d:\kmer.db d:\matrix.csv
one2all d:\kmer.db E:\Data\kmc250\GCF_000171975.1_ASM17197v1_genomic_s d:\vector.csv
distance d:\matrix.csv
list-patterns d:\kmer.db d:\patterns.txt

minhash E:\Data\kmc250.list 0.1
build-mh E:\Data\kmc250.list d:\kmer-mh.db
all2all d:\kmer-mh.db d:\matrix-mh.csv
one2all d:\kmer-mh.db E:\Data\kmc250\GCF_000171975.1_ASM17197v1_genomic_s d:\vector-mh.csv



*/


int main(int argc, char **argv)
{
	Console console;
	console.parse(argc, argv);
	


	return 0;
}

