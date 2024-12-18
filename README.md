# Kmer-db
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/kmer-db/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/kmer-db/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/kmer-db.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/kmer-db)
[![Build and tests](../../workflows/Build%20and%20tests/badge.svg)](../../actions/workflows/main.yml)
[![License](https://anaconda.org/bioconda/famsa/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M&color=yellow&logo=Apple&logoColor=white)
![Windows](https://img.shields.io/badge/%E2%80%8B-Windows-00A98F?logo=windows)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)

Kmer-db is a fast and memory-efficient tool for large-scale k-mer analyses (indexing, querying, estimating evolutionary relationships, etc.). 

## Quick start

```bash
git clone --recurse-submodules https://github.com/refresh-bio/kmer-db
cd kmer-db && gmake

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

Kmer-db comes with a set of [precompiled binaries](https://github.com/refresh-bio/kmer-db/releases) for Linux, macOS, and Windows. 
The software is also available on [Bioconda](https://anaconda.org/bioconda/kmer-db):
```
conda install -c bioconda kmer-db
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/).
Kmer-db can be also built from the sources distributed as:

* GNU Make project for Linux and macOS (gmake 4.3 and gcc/g++ 11 or newer required),
* Visual Studio 2022 solution for Windows.


## Vector extensions

Kmer-db can be built for x86-64 and ARM64 8 architectures (including Apple Mx based on ARM64 8.4 core) and takes advantage of AVX2 (x86-64) and NEON (ARM) CPU extensions. The default target platform is x86-64 with AVX2 extensions. This, however, can be changed by setting `PLATFORM` variable for `make`:

```bash
make PLATFORM=none    # unspecified platform, no extensions
make PLATFORM=sse2    # x86-64 with SSE2 
make PLATFORM=avx     # x86-64 with AVX 
make PLATFORM=avx2    # x86-64 with AVX2 (default)
make PLATFORM=native  # x86-64 with AVX2 and native architecture
make PLATFORM=arm8    # ARM64 8 with NEON  
make PLATFORM=m1      # ARM64 8.4 (especially Apple M1) with NEON 
```   

Note, that x86-64 binaries determine the supported extensions at runtime, which makes them backwards-compatible. For instance, the AVX executable will also work on SSE-only platform, but with limited performance.

# 2. Usage

`kmer-db <mode> [options] <positional arguments>`

Kmer-db operates in one of the following modes:

* `build` - building a database from samples,
* `all2all` - counting common k-mers - all samples in the database,
* `all2all-sp` - counting common k-mers - all samples in the database (sparse computation),
* `all2all-parts` - counting common k-mers - all samples within from databases (sparse computation),
* `new2all` - counting common k-mers - set of new samples versus database,
* `one2all` - counting common k-mers - single sample versus database,
* `distance` - calculating similarities/distances,
* `minhash` - storing minhashed k-mers.
    
Common options:
* `-t <threads>` - number of threads (default: number of available cores),
    
The meaning of other options and positional arguments depends on the selected mode.
    
## 2.1. Building a database

Construction of k-mers database is an obligatory step for further analyses. The procedure accepts several input types:
* compressed or uncompressed genomes/reads:

    ```kmer-db build [-k <kmer-length>] [-f <fraction>] [-multisample-fasta] [-extend] [-alphabet <alphabet>] [-preserve-strand] [-t <threads>] <samples> <database>```
* [KMC-generated](https://github.com/refresh-bio/KMC) k-mers: 

    ```kmer-db build -from-kmers [-f <fraction>] [-extend] [-t <threads>] <samples> <database>```
  
* [minhashed k-mers](#24-storing-minhashed-k-mers) produced by `minhash` mode:

    ```kmer-db build -from-minhash [-extend] [-t <threads>] <samples> <database>```

Parameters:
 * `samples` (input) - one of the following: 
    * FASTA file (*fa*, *fna*, *fasta*, *fa.gz*, *fna.gz*, *fasta.gz*) with one or multiple (`-multisample-fasta` switch) samples
    * file with a newline-separated list of samples:
      ```
      sample_file_1
      sample_file_2
      sample_file_3
      ...
      ```
   Every file can be in one of the formats:
     1. FASTA genomes/reads (default). If a file on the list cannot be found, the following extensions are tested: *fa*, *fna*, *fasta*, *gz*, *fa.gz*, *fna.gz*, *fasta.gz*.
     2. [KMC-generated](https://github.com/refresh-bio/KMC) k-mer files (`-from-kmers` switch specified). A set of two KMC files (*.kmc_pre* + *.kmc_suf*) is required for every list entry.
     3. minhashed k-mers (`-from-minhash` switch specified). Minhashed k-mer files (*.minhash*) must be generated by `minhash` command [prior to the database construction](#24-storing-minhashed-k-mers).  
	Note, that minhashing may be also done during the database construction by specyfying `-f` option.
* `database` (output) - file with generated k-mer database, 
* `-k <kmer-length>` - length of k-mers (default: 18); ignored when `-from-kmers` or `-from-minhash` switch is specified,
* `-f <fraction>` - fraction of all k-mers to be accepted by the minhash filter during database construction (default: 1); ignored when `-from-minhash` switch is present,
* `-multisample-fasta` - each sequence in a FASTA file is treated as a separate sample,
* `-extend` - extend the existing database with new samples,
* `-alphabet` - alphabet:
  * `nt` (4 symbol nucleotide with indistinguishable T/U; default) 
  * `aa` (20 symbol amino acid)
  * `aa12_mmseqs` (amino acid reduced to 12 symbols as in MMseqs: AST,C,DN,EQ,FY,G,H,IV,KR,LM,P,W
  * `aa11_diamond` (amino acid reduced to 11 symbols as in Diamond: KREDQN,C,G,H,ILV,M,F,Y,W,P,STA
  * `aa6_dayhoff` (amino acid reduced to 6 symbols as proposed by Dayhoff: STPAG,NDEQ,HRK,MILV,FYW,C
* `-preserve-strand`- preserve strand instead of taking canonical k-mers (allowed only in `nt` alphabet; default: off)
* `-t <threads>` - number of threads (default: number of available cores).

## 2.2. Counting common k-mers 

### Samples in the database against each other:

Dense computations - recomended when the distance matrix contains few zeros. Output can be stored in the dense or sparse form (`-sparse` switch).

`kmer-db all2all [-buffer <size_mb>] [-t <threads>] [-sparse [-min [<criterion>:]<value>]* [-max [<criterion>:]<value>]* ] <database> <common_table>`
 
Sparse computations - recommended when the distance matrix contains many zeros. Output matrix is always in the sparse form:

`kmer-db all2all-sp [-buffer <size_mb>] [-t <threads>] [-min [<criterion>:]<value>]* [-max [<criterion>:]<value>]* [-sample-rows [<criterion>:]<count>] <database> <common_table>`

Sparse computations, partial databases - use when the distance matrix contains many zeros and there are multiple partial databases. Output matrix is always in the sparse form:

`kmer-db all2all-parts [-buffer <size_mb>] [-t <threads>] [-min [<criterion>:]<value>]* [-max [<criterion>:]<value>]* [-sample-rows [<criterion>:]<count>] <db_list> <common_table>`
 
Parameters:
* `database` (input) - k-mer database file created by `build` mode,
* `db_list` (input) - file containing list of databases files created by `build` mode,
* `common_table` (output) - file containing table with common k-mer counts,
* `-buffer <size_mb>` - size of cache buffer in megabytes; use L3 size for Intel CPUs and L2 for AMD for best performance; default: 8,
* `-t <threads>` - number of threads (default: number of available cores),
* `-sparse` - stores output matrix in a sparse form (always on in `all2all-sp` and `all2all-parts` modes),
* `-min [<criterion>:]<value>` - retains elements with `criterion` greater than or equal to `value` (see details below), 
* `-max [<criterion>:]<value>` - retains elements with `criterion` lower than or equal to `value` (see details below),
* `-sample-rows [<criterion>:]<count>` - retains `count` elements in every row using one of the strategies: (i) random selection (no `criterion`); (ii) the best elements with respect to `criterion`.

`criterion` can be `num-kmers` (number of common k-mers) or one of the distance/similarity measures: `jaccard`, `min`, `max`, `cosine`, `mash`, `ani`, `ani-shorder` (see 2.3 for definitions). No `criterion` indicates `num-kmers` (filtering) or random elements selection (sampling). Multiple filters can be combined. 

### New samples against the database:

`kmer-db new2all [-multisample-fasta | -from-kmers | -from-minhash] [-t <threads>]  [-sparse [-min [<criterion>:]<value>]* [-max [<criterion>:]<value>]* ] <database> <samples> <common_table>`

Parameters:
* `database` (input) - k-mer database file created by `build` mode,
* `samples` (input) - file containing samples in one of the supported formats (see `build` mode); if samples are given as genomes (default) or k-mers (`-from-kmers` switch), the minhashing is done automatically with the same filter as in the database,
* `common_table` (output) - file containing table with common k-mer counts,
* `-multisample-fasta` / `-from-kmers` / `-from-minhash` - see `build` mode for details,
* `-t <threads>` - number of threads (default: number of available cores),
* `-sparse` - stores output matrix in a sparse form,
* `-min [<criterion>:]<value>` - retains elements with `criterion` greater than or equal to `value` (see details below), 
* `-max [<criterion>:]<value>` - retains elements with `criterion` lower than or equal to `value` (see details below),

`criterion` can be `num-kmers` (number of common k-mers) or one of the distance/similarity measures: `jaccard`, `min`, `max`, `cosine`, `mash`, `ani`, `ani-shorder` (see 2.3 for definitions). No `criterion` indicates `num-kmers`. Multiple filters can be combined.
 
### Single sample against the database:

`kmer-db one2all [-from-kmers | -from-minhash] [-t <threads>] <database> <sample> <common_table>`

The meaning of the parameters is the same as in `new2all` mode, but instead of specifying file with sample list, a single sample file is used as a query.

### Output format

Modes `all2all`, `all2all-sp`, `all2all-parts`, `new2all`, and `one2all` produce a comma-separated table with numbers of common k-mers. For `all2all`, `new2all`, and `one2all` modes, the table is by default stored in a dense form:

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

When `-sparse` switch is specified or `all2all-sp`, `all2all-parts` modes are used, the table is stored in a sparse form. In particular, zeros are omitted while non-zero elements are represented as pairs (*column_id*: *value*) with 1-based column indexing. Thus, rows may have different number of elements, e.g.:

| 									| 								| 					| 				|		|			|	
| :---: 							| :---: 						| :---: 			| :---:			| :---:	|  :---:	| 
| kmer-length: *k* fraction: *f* 	| db-samples 					| *s<sub>1</sub>*					| *s<sub>2</sub>* | ... 	|  *s<sub>n</sub>* |
| query-samples 					| total-kmers 					| &#124;*s<sub>1</sub>*&#124;		| &#124;*s<sub>2</sub>*&#124; 	| ... 	|  &#124;*s<sub>n</sub>*&#124; |
| *q<sub>1</sub>* 					| &#124;*q<sub>1</sub>*&#124;	| *i<sub>11</sub>*: &#124;*q<sub>1</sub> &cap; s<sub>i<sub>11</sub></sub>*&#124;	| *i<sub>12</sub>*: &#124;*q<sub>1</sub> &cap; s<sub>i<sub>12</sub></sub>*&#124; | ||
| *q<sub>2</sub>* 					| &#124;*q<sub>2</sub>*&#124;	| *i<sub>21</sub>*: &#124;*q<sub>2</sub> &cap; s<sub>i<sub>21</sub></sub>*&#124;	| *i<sub>22</sub>*: &#124;*q<sub>2</sub> &cap; s<sub>i<sub>22</sub></sub>*&#124; 	| *i<sub>23</sub>*: &#124;*q<sub>2</sub> &cap; s<sub>i<sub>23</sub></sub>*&#124;  	| |   
| *q<sub>2</sub>* 					| &#124;*q<sub>2</sub>*&#124;	| ||||
| ... 								| ...							| ... ||||
| *q<sub>m</sub>* 					| &#124;*q<sub>m</sub>*&#124;	| *i<sub>m1</sub>*: &#124;*q<sub>m</sub> &cap; s<sub>i<sub>m1</sub></sub>*&#124;	| |||

For performance reasons, `all2all`, `all2all-sp`, and `all2all-parts` modes produce a lower triangular matrix.

 
 ## 2.3. Calculating similarities or distances

`kmer-db distance <measure> [-sparse [-min [<criterion>:]<value>]* [-max [<criterion>:]<value>]* ] <common_table> <output_table>`

Parameters:
* `measure` - names of the similarity/distance measure to be calculated, can be one of the following: 
  * `jaccard`: $J(q,s) = |p \cap q| / |p \cup q|$, 
  * `min`: $\min(q,s) =  |p \cap q| / \min(|p|,|q|)$, 
  * `max`: $\max(q,s) =  |p \cap q| / \max(|p|,|q|)$, 
  * `cosine`: $\cos(q,s) = |p \cap q| / \sqrt{|p| \cdot |q|}$,  
  * `mash` (Mash distance): $\textrm{Mash}(q,s) = -\frac{1}{k}ln\frac{2 \cdot J(q,s)}{1 + J(q,s)}$, 
  * `ani` (average nucleotide identity): $\textrm{ANI}(q,s) = 1 - \textrm{Mash}(p,q)$,
  * `ani-shorter` - same as `ani` but with `min` used instead of `jaccard`.
* `common_table` (input) - file containing table with numbers of common k-mers produced by `all2all`, `new2all`, or `one2all` mode (both, dense and sparse matrices are supported), 
* `output_table` (output) - file containing table with calculated distance measure,  
* `-phylip-out` - store output distance matrix in a Phylip format,
* `-sparse` - outputs a sparse matrix (only for dense input matrices - sparse inputs always produce sparse outputs),
* `-min [<criterion>:]<value>` - retains elements with `criterion` greater than or equal to `value` (see details below), 
* `-max [<criterion>:]<value>` - retains elements with `criterion` lower than or equal to `value` (see details below),

`criterion` can be `num-kmers` (number of common k-mers) or one of the distance/similarity measures: `jaccard`, `min`, `max`, `cosine`, `mash`, `ani`, `ani-shorder` (see 2.3 for definitions). If no `criterion` is specified, `measure` argument is used by default. Multiple filters can be combined.
 
        
## 2.4. Storing minhashed k-mers

This is an optional analysis step which stores minhashed k-mers on the hard disk to be later consumed by `build`, `new2all`, or `one2all` modes with `-from-minhash` switch. It can be skipped if one wants to use all k-mers from samples for distance estimation or employs minhashing during database construction. Syntax:

`kmer-db minhash [-f <fraction>] [-k <kmer-length>] [-multisample-fasta] [-alphabet <alphabet>] [-preserve-strand] <samples>`

`kmer-db minhash -from-kmers [-f <fraction>] <samples>`

Parameters:
 * `sample_list` (input) - file containing list of samples in one of the supported formats (see `build` mode), 
 * `-f <fraction>` - fraction of all k-mers to be accepted by the minhash filter (default: 0.01),
 * `-k <kmer-length>` - length of k-mers (default: 18; maximum: 30); ignored when `-from-kmers` switch is specified,
 * `-multisample-fasta` / `-from-kmers` - see `build` mode for details.
 * `-alphabet` - alphabet:
   * `nt` (4 symbol nucleotide with indistinguishable T/U; default) 
   * `aa` (20 symbol amino acid)
   * `aa12_mmseqs` (amino acid reduced to 12 symbols as in MMseqs: AST,C,DN,EQ,FY,G,H,IV,KR,LM,P,W
   * `aa11_diamond` (amino acid reduced to 11 symbols as in Diamond: KREDQN,C,G,H,ILV,M,F,Y,W,P,STA
   * `aa6_dayhoff` (amino acid reduced to 6 symbols as proposed by Dayhoff: STPAG,NDEQ,HRK,MILV,FYW,C
* `-preserve-strand`- preserve strand instead of taking canonical k-mers (allowed only in `nt` alphabet; default: off)

For each sample from the list, a binary file with *.minhash* extension containing filtered k-mers is created.


# 3. Datasets
List of the pathogens investigated in Kmer-db study can be found [here](https://github.com/refresh-bio/kmer-db/tree/master/data)

## Citing
[Deorowicz, S., Gudyś, A., Długosz, M., Kokot, M., Danek, A. (2019) Kmer-db: instant evolutionary distance estimation, Bioinformatics, 35(1): 133–136](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty610/5050791?redirectedFrom=fulltext)
