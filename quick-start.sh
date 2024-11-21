# !/bin/bash
INPUT=./test/virus
OUTPUT=./output

mkdir $OUTPUT

# build a database from all 18-mers (default) contained in a set of sequences
./bin/kmer-db build $INPUT/seqs.part1.list $OUTPUT/k18.db

# establish numbers of common k-mers between new sequences and the database
./bin/kmer-db new2all $OUTPUT/k18.db $INPUT/seqs.part2.list $OUTPUT/n2a.csv

# calculate jaccard index from common k-mers
./bin/kmer-db distance jaccard $OUTPUT/n2a.csv $OUTPUT/n2a.jaccard

# extend the database with new sequences
./bin/kmer-db build -extend $INPUT/seqs.part2.list $OUTPUT/k18.db

# establish numbers of common k-mers between all sequences in the database
./bin/kmer-db all2all $OUTPUT/k18.db $OUTPUT/a2a.csv

# build a database from 10% of 25-mers using 16 threads
./bin/kmer-db build -k 25 -f 0.1 -t 16 $INPUT/seqs.part1.list $OUTPUT/k25.db

# establish number of common 25-mers between single sequence and the database 
# (minhash filtering that retains 10% of MT159713 k-mers is done prior to the comparison)  
./bin/kmer-db one2all $OUTPUT/k25.db $INPUT/data/MT159713.fasta $OUTPUT/MT159713.csv

# build two partial databases
./bin/kmer-db build $INPUT/seqs.part1.list  $OUTPUT/k18.parts1.db
./bin/kmer-db build $INPUT/seqs.part2.list  $OUTPUT/k18.parts2.db

# establish numbers of common k-mers between all sequences in the databases,
# computations are done in the sparse mode, the output matrix is also sparse
echo $OUTPUT/k18.parts1.db > $OUTPUT/db.list
echo $OUTPUT/k18.parts2.db >> $OUTPUT/db.list
./bin/kmer-db all2all-parts $OUTPUT/db.list $OUTPUT/k18.parts.csv
