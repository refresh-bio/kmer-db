# count_kmc_databases.py

count_kmc_databases.py is a tool which performs *k*-mer counting for many small FASTA files.
It simultaneously runs many single-threaded KMC instances to speedup the counting.

## Usage
Firstly, it is necessary to define the following internal script variables:
* `kmc` - path to KMC binary,
* `kmc_tools` - path to KMC tools binary.

Next, it can by run with the following command:

`count_kmc_databases.py <k> <threads> <list_of_fasta>`

Parameters:
* `k` - *k*-mer length,
* `threads` - number of threads,
* `list_of_fasta` - file containing list of genomes files in .fna.gz format.
