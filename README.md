# kmer-db-dev

## USAGE
`kmer-db <mode> [options] <positional arguments>`

KmerDB operates in one of the following modes:

* `minhash` - minhashing k-mers,
* `build` - building a database from k-mers,
* `all2all` - calculating number of common k-mers between all samples in the database,
* `one2all` - calculating number of common kmers between single sample and all the samples in the database,
* `distance` - calculating similarities/distances.
    
Options:

* `-t num-threads` - distributes processing over `num-threads`,
* `-mh-input` - loads samples from minhashed files (applies only to `build` and `one2all` modes).
    
The meaning of the positional arguments depends on the selected mode.
    
### Minhashing k-mers
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
    When reading k-mers, the tool searches for KMC files with `.pre` and `.suf` extensions for each sample. 
 * `fraction` (input) - a number from <0,1> interval determining a fraction of all k-mers to be accepted by the filter.
 
  For each sample from the list, a binary file with `.minhash` extension containing filtered k-mers is created.

### Building the database from k-mers
Construction of k-mers database is an obligatory step for further analyses. Both, raw [KMC-generated](https://github.com/refresh-bio/KMC), as well as minhashed k-mers are accepted as an input for the procedure.  

```
kmer-db build <sample_list> <database>
kmer-db build -mh-input <sample_list> <database>
```

Parameters:
 * `sample_list` (input) - file containing list of k-mer files, one entry per sample. The tool requires corresponding KMC-generated files (with `.pre` and `.suf` extensions), or minhashed files (`.minhash`) if `-mh-input` switch is specified.
 * `database` (output) - file with k-mer database for all samples.
 
 ### Calculating number of common k-mers ###
Calculating number of common k-mers for all the samples in the database:
 
 `kmer_db all2all <database> <common_matrix>`
 
Parameters:
* `database` (input) - k-mer database file created by `build` mode,
* `common_matrix` (output) - file containing matrix with common k-mer counts.

Calculating number of common kmers between single sample and all the samples in the database:

`kmer_db one2all <database> <sample> <common_vector>`
`kmer_db one2all -mh-input <database> <sample> <common_vector>`

Parameters:
 * `database` (input) - k-mer database file created by `build` mode,
 * `sample` (input) - k-mer file for a sample (raw or minhashed, depending on the presence of `-mh-input` switch),
 * `common_vector` (output) - file containing vector with numbers of common k-mers.
 
 ### Calculating similarities/distances
 
`kmer_db distance <common_file>"`

Parameters:
* `common_file` (input) - file containing matrix (vector) with numbers of common k-mers produced by `all2all` (`one2all`) mode.

This mode generates a set of files containing matrices (vectors) with different similarity/distance measures:
* `<common_file>.jaccard`
* `<common_file>.min` 
* `<common_file>.max` 
* `<common_file>.cosine` 
* `<common_file>.mash`


