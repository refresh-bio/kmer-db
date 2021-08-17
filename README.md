# Kmer-db
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/kmer-db/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/kmer-db/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/kmer-db.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/kmer-db)
[![C/C++ CI](https://github.com/refresh-bio/kmer-db/workflows/C/C++%20CI/badge.svg)](https://github.com/refresh-bio/kmer-db/actions)


Kmer-db is a fast and memory-efficient tool for large-scale k-mer analyses (indexing, querying, estimating evolutionary relationships, etc.). 

## Quick start

```bash
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
```


### Table of contents

1. [Installation](#1-installation)
2. [Usage](#2-usage)
    1. [Building a database](#21-building-a-database)
    2. [Counting common k-mers](#22-counting-common-k-mers)
    3. [Calculating similarities or distances](#23-calculating-similarities-or-distances)
    4. [Storing minhashed k-mers](#24-storing-minhashed-k-mers)
3. [Datasets](#3-datasets)


# 1. Installation

Kmer-db comes with a set of [precompiled binaries](https://github.com/refresh-bio/kmer-db/releases) for Linux, OS X, and Windows. 
The software is also available on [Bioconda](https://anaconda.org/bioconda/kmer-db):
```
conda install -c bioconda kmer-db
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda).
Kmer-db can be also built from the sources distributed as:

* MAKE project (G++ 4.8.5 tested) for Linux and OS X,
* Visual Studio 2015 solution for Windows.



## *zlib* linking

Kmer-db uses *zlib* for handling gzipped inputs. Under Linux, the software is by default linked against system-installed *zlib*. Due to issues with some library versions, precompiled *zlib* is also present the repository. In order to use it, one needs to modify variable INTERNAL_ZLIB at the top of the makefile. Under Windows, the repository library is always used.

## AVX and AVX2 support

Kmer-db, by default, takes advantage of AVX (required) and AVX2 (optional) CPU extensions. The pre-built binary determines supported instructions at runtime, thus it is multiplatform. When compiling the sources under Linux and OS X, the support of AVX2 is also established automatically. Under Windows, the program is by default built with AVX2 instructions. To prevent this, Kmer-db must be compiled with NO_AVX2 symbolic constant defined.

# 2. Usage

`kmer-db <mode> [options] <positional arguments>`

Kmer-db operates in one of the following modes:

* `build` - building a database from samples,
* `all2all` - counting common k-mers - all samples in the database,
* `new2all` - counting common k-mers - set of new samples versus database,
* `one2all` - counting common k-mers - single sample versus database,
* `distance` - calculating similarities/distances,
* `minhash` - storing minhashed k-mers,
    
Common options:
* `-t <threads>` - number of threads (default: number of available cores),
    
The meaning of other options and positional arguments depends on the selected mode.
    
## 2.1. Building a database

Construction of k-mers database is an obligatory step for further analyses. The procedure accepts several input types:
* compressed or uncompressed genomes:

    ```kmer-db build [-k <kmer-length>] [-f <fraction>] [-multisample-fasta] [-extend] <sample_list> <database>```
* [KMC-generated](https://github.com/refresh-bio/KMC) k-mers: 

    ```kmer-db build -from-kmers [-f <fraction>] [-extend] <sample_list> <database>```
  
* [minhashed k-mers](#24-storing-minhashed-k-mers) produced by `minhash` mode:

    ```kmer-db build -from-minhash [-extend] <sample_list> <database>```

Parameters:
* `sample_list` (input) - file containing list of samples in the following format:
    ```
    sample_file_1
    sample_file_2
    sample_file_3
    ...
    ```
    By default, the tool requires uncompressed or compressed FASTA files for each sample. If a file on the list cannot be found, the package tries adding the following extensions: *fna*, *fasta*, *gz*, *fna.gz*, *fasta.gz* . When `-from-kmers` switch is specified, corresponding [KMC-generated](https://github.com/refresh-bio/KMC) k-mer files (*.kmc_pre* and *.kmc_suf*) are required. If `-from-minhash` switch is present, minhashed k-mer files (*.minhash*) must be generated by `minhash` command [prior to the database construction](#24-storing-minhashed-k-mers). Note, that minhashing may be also done during the database construction by specyfying `-f` option.
* `database` (output) - file with generated k-mer database. 
* `-k <kmer-length>` - length of k-mers (default: 18); ignored when `-from-kmers` or `-from-minhash` switch is specified.
* `-f <fraction>` - fraction of all k-mers to be accepted by the minhash filter during database construction (default: 1); ignored when `-from-minhash` switch is present.
* `-multisample-fasta` - each sequence in a genome FASTA file is treated as a separate sample,
* `-extend` - extend the existing database with new samples.

## 2.2. Counting common k-mers 

### Samples in the database against each other:
 
 `kmer-db all2all [-buffer <size_mb>] [-sparse] <database> <common_table>`
 
Parameters:
* `database` (input) - k-mer database file created by `build` mode,
* `common_table` (output) - file containing table with common k-mer counts.
* `-buffer <size_mb>` - size of cache buffer in megabytes; use L3 size for Intel CPUs and L2 for AMD for best performance; default: 8
* `-sparse` - stores output matrix in a sparse form.

### New samples against the database:

`kmer-db new2all [-multisample-fasta | -from-kmers | -from-minhash] [-sparse] <database> <sample_list> <common_table>`

Parameters:
* `database` (input) - k-mer database file created by `build` mode.
* `sample_list` (input) - file containing list of samples in one of the supported formats (see `build` mode); if samples are given as genomes (default) or k-mers (`-from-kmers` switch), the minhashing is done automatically with the same filter as in the database.
* `common_table` (output) - file containing table with common k-mer counts.
* `-multisample-fasta` / `-from-kmers` / `-from-minhash` - see `build` mode for details.
* `-sparse` - stores output matrix in a sparse form.
 
### Single sample against the database:

`kmer-db one2all [-multisample-fasta | -from-kmers | -from-minhash] <database> <sample> <common_table>`

The meaning of the parameters is the same as in `new2all` mode, but instead of specifying file with sample list, a single sample file is used as a query.

### Output format

Modes `all2all`, `new2all`, and `one2all` produce a comma-separated table with numbers of common k-mers. The table is by default stored in a dense form:

| 									| 								| 					| 					|		|			|	
| :---: 							| :---: 						| :---: 			| :---:				| :---:	|  :---:	| 
| kmer-length: *k* fraction: *f* 	| db-samples 					| *s<sub>1</sub>*									| *s<sub>2</sub>* 										| ... 	|  *s<sub>n</sub>* |
| query-samples 					| total-kmers 					| &#124;*s<sub>1</sub>*&#124;						| &#124;*s<sub>2</sub>*&#124; 							| ... 	|  &#124;*s<sub>n</sub>*&#124; |
| *q<sub>1</sub>* 					| &#124;*q<sub>1</sub>*&#124;	| &#124;*q<sub>1</sub> &cap; s<sub>1</sub>*&#124;	| &#124;*q<sub>1</sub> &cap; s<sub>2</sub>*&#124; 	| ... 	|  &#124;*q<sub>1</sub> &cap; s<sub>n</sub>*&#124; |
| *q<sub>2</sub>* 					| &#124;*q<sub>2</sub>*&#124;	| &#124;*q<sub>2</sub> &cap; s<sub>1</sub>*&#124;	| &#124;*q<sub>2</sub> &cap; s<sub>2</sub>*&#124; 	| ... 	|  &#124;*q<sub>2</sub> &cap; s<sub>n</sub>*&#124; |
| ... 								| ...							| ...												| ...													| ...	|	...												 |
| *q<sub>m</sub>* 					| &#124;*q<sub>m</sub>*&#124;	| &#124;*q<sub>m</sub> &cap; s<sub>1</sub>*&#124;	| &#124;*q<sub>m</sub> &cap; s<sub>2</sub>*&#124; 	| ... 	|  &#124;*q<sub>m</sub> &cap; s<sub>n</sub>*&#124; |

where:
* *k* - k-mer length,
* *f* - minhash fraction (1, when minhashing is disabled), 
* *s<sub>1</sub>*, *s<sub>2</sub>*,  ...,   *s<sub>n</sub>* - database sample names,
* *q<sub>1</sub>*, *q<sub>2</sub>*,  ...,   *q<sub>m</sub>* - query sample names,
* &#124;*a*&#124; - number of k-mers in sample *a*,
* &#124;*a &cap; b*&#124; - number of k-mers common for samples *a* and *b*.

For performance reasons, `all2all` mode produces a lower triangular matrix.

When `-sparse` switch is specified, the table is stored in a sparse form. In particular, zeros are omitted while non-zero elements are represented as pairs (*column_id*, *value*) with 1-based column indexing. Thus, rows may have different number of elements, e.g.:

| 									| 								| 					| 				|		|			|	
| :---: 							| :---: 						| :---: 			| :---:			| :---:	|  :---:	| 
| kmer-length: *k* fraction: *f* 	| db-samples 					| *s<sub>1</sub>*					| *s<sub>2</sub>* | ... 	|  *s<sub>n</sub>* |
| query-samples 					| total-kmers 					| &#124;*s<sub>1</sub>*&#124;		| &#124;*s<sub>2</sub>*&#124; 	| ... 	|  &#124;*s<sub>n</sub>*&#124; |
| *q<sub>1</sub>* 					| &#124;*q<sub>1</sub>*&#124;	| (*i<sub>11</sub>*, &#124;*q<sub>1</sub> &cap; s<sub>i<sub>11</sub></sub>*&#124;)	| (*i<sub>12</sub>*, &#124;*q<sub>1</sub> &cap; s<sub>i<sub>12</sub></sub>*&#124;) | ||
| *q<sub>2</sub>* 					| &#124;*q<sub>2</sub>*&#124;	| (*i<sub>21</sub>*, &#124;*q<sub>2</sub> &cap; s<sub>i<sub>21</sub></sub>*&#124;)	| (*i<sub>22</sub>*, &#124;*q<sub>2</sub> &cap; s<sub>i<sub>22</sub></sub>*&#124;) 	| (*i<sub>23</sub>*, &#124;*q<sub>2</sub> &cap; s<sub>i<sub>23</sub></sub>*&#124;)  	| |   
| *q<sub>2</sub>* 					| &#124;*q<sub>2</sub>*&#124;	| ||||
| ... 								| ...							| ... ||||
| *q<sub>m</sub>* 					| &#124;*q<sub>m</sub>*&#124;	| (*i<sub>m1</sub>*, &#124;*q<sub>m</sub> &cap; s<sub>i<sub>m1</sub></sub>*&#124;)	| |||
 
 ## 2.3. Calculating similarities or distances

`kmer-db distance [<measures>] <common_table>`

Parameters:
* `common_table` (input) - file containing table with numbers of common k-mers produced by `all2all`, `new2all`, or `one2all` mode (only dense matrices are supported). 
* `measures` - names of the similarity/distance measures to be calculated, can be one or several of the following: `jaccard`, `min`, `max`, `cosine`, `mash`. If measures are not specified, `jaccard` is used by default.
* `-phylip-out` - store output distance matrix in a Phylip format.

This mode generates a file with similarity/distance table for each selected measure. Name of the output file is produced by adding to the input file an extension with a measure name.
    
    
## 2.4. Storing minhashed k-mers

This is an optional analysis step which stores minhashed k-mers on the hard disk to be later consumed by `build`, `new2all`, or `one2all` modes with `-from-minhash` switch. It can be skipped if one wants to use all k-mers from samples for distance estimation or employs minhashing during database construction. Syntax:

`kmer-db minhash [-k <kmer-length>] [-multisample-fasta] <fraction> <sample_list>`

`kmer-db minhash -from-kmers <fraction> <sample_list>`

Parameters:
 * `fraction` (input) - fraction of all k-mers to be accepted by the minhash filter.
 * `sample_list` (input) - file containing list of samples in one of the supported formats (see `build` mode). 
 * `-k <kmer-length>` - length of k-mers (default: 18); ignored when `-from-kmers` switch is specified.
 * `-multisample-fasta` / `-from-kmers` - see `build` mode for details.

For each sample from the list, a binary file with *.minhash* extension containing filtered k-mers is created.


# 3. Datasets
List of the pathogens investigated in Kmer-db study can be found [here](https://github.com/refresh-bio/kmer-db/tree/master/data)

## Citing
[Deorowicz, S., Gudyś, A., Długosz, M., Kokot, M., Danek, A. (2019) Kmer-db: instant evolutionary distance estimation, Bioinformatics, 35(1): 133–136](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty610/5050791?redirectedFrom=fulltext)
