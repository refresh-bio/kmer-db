# !/bin/bash
git clone https://github.com/refresh-bio/kmer-db
cd kmer-db && make

INPUT=./test/virus
OUTPUT=./output

mkdir $OUTPUT

# build a database from all 18-mers (default) contained in a set of sequences
./kmer-db build $INPUT/seqs.part1.list $OUTPUT/k18.db

# establish numbers of common k-mers between new sequences and the database
./kmer-db new2all $OUTPUT/k18.db $INPUT/seqs.part2.list $OUTPUT/n2a.csv

# calculate jaccard index from common k-mers
./kmer-db distance $OUTPUT/n2a.csv

# extend the database with new sequences
./kmer-db build -extend $INPUT/seqs.part2.list $OUTPUT/k18.db

# establish numbers of common k-mers between all sequences in the database
./kmer-db all2all $OUTPUT/k18.db $OUTPUT/a2a.csv

# build a database from 10% of 25-mers using 16 threads
./kmer-db build -k 25 -f 0.1 -t 16 $INPUT/seqs.part1.list $OUTPUT/k25.db

# establish number of common 25-mers between single sequence and the database 
# (minhash filtering that retains 10% of MT159713 k-mers is done prior to the comparison)  
./kmer-db one2all $OUTPUT/k25.db $INPUT/data/MT159713.fasta $OUTPUT/MT159713.csv