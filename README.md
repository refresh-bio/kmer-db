# kmer-db-dev

## USAGE

### Min-hashing k-mers
This is an optional analysis step - it can be skipped if one wants to use all k-mers from samples for distance estimation. 

`kmer-db minhash <sample_list> <fraction>`

Parameters:
 * `sample_list` (input) - file containing list of [KMC-generated](https://github.com/refresh-bio/KMC) k-mer files for samples to be analyzed. In the list file there should be a single entry for each sample:
    ```
    sample1
    sample2
    sample3
    ...
    ```
    When reading k-mers, the tool will automatically add `.pre` and `.suf` extensions for each sample. 
 * `fraction` (input) - a number from <0,1> interval determining a fraction of all k-mers to be accepted by filter.
 
  For each sample from the list, a binary file with `.minhash` extension containing filtered k-mers is created.

### Building the database from k-mers
Construction of k-mers database is an obligatory step for further analyses. Both, raw [KMC-generated](https://github.com/refresh-bio/KMC) as well as min-hashed k-mers are accepted as in input for the procedure.  

`kmer-db build <sample_list> <database>`
`kmer-db build-minhash <sample_list> <database>`

Parameters:
 * `sample_list` (input) - file containing list of k-mer files, one entry per sample. In `build` mode, the tool requires corresponding KMC-generated files (with `.pre` and `.suf` extensions), while in `build-minhash`, min-hashed files (`.minhash`) are neccessary.
 * `database` (output) - file with k-mer database for all samples.
 
 ### Calculating number of common k-mers ###
Calculating number of common k-mers for all the samples in a database:
 
 `kmer_db all2all <database> <common_matrix>`
 
Parameters:
* `database` (input) - k-mer database file,
* `common_matrix` (output) - file containing matrix with common k-mer counts

Calculating number of common kmers between single sample and all the samples in the database:

`kmer_db one2all <database> <sample> <common_vector>`

Parameters:
 * `database` (input) - k-mer database file
 * `sample` (input) - k-mer file for a sample (raw or min-hashed)
 * `common_vector` (output) - file containing vector with numbers of common k-mers
 
 ### Calculating distance metrices/vectors
 
`kmer_db distance <common_file>"`

Parameters:
* `common_file` (input) - file containing matrix/vector with numbers of common k-mers

This mode generates set of files with matrices/vectors with different similarity measures:
* `<common_file>.jaccard` 
* `<common_file>.mash`

