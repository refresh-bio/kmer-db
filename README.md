# kmer-db-dev

## USAGE
`kmer-db <mode> [options] <positional arguments>`

Kmer-db operates in one of the following modes:

* `minhash` - minhashing k-mers,
* `build` - building a database from k-mers,
* `build-mh` - building a database from minhashed k-mers,
* `all2all` - calculating number of common k-mers between all samples in the database,
* `one2all` - calculating number of common kmers between single sample and all the samples in the database,
* `distance` - calculating similarities/distances.
    
Options:

* `-t <threads>` - number of threads (default: number of available cores),
* `-c <size_mb>` - size of cache buffer in megabytes, applies to `all2all` mode (default: 8; use L2 size for Intel CPUs and L3 for AMD to maximize performance).
    
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
kmer-db build-mh <sample_list> <database>
```

Parameters:
 * `sample_list` (input) - file containing list of k-mer files, one entry per sample. In `build` mode, the tool requires corresponding KMC-generated files (with `.pre` and `.suf` extensions). In `build-mh`, minhashed files are required (`.minhash`).
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

### Toy examples

Let `pathogens.list` be the file containing names of KMC-processed samples (there exist `.pre` and `.suf` files for each sample):
```
acinetobacter
klebsiella
e.coli
...
```

Calculating similarities/distances between all samples listed in `pathogens.list` using all k-mers. 
```
kmer-db build pathogens.list pathogens.db
kmer-db all2all pathogens.db matrix.csv
kmer-db distance matrix.csv
```

Same as above, but using only 10% of k-mers.
```
kmer-db minhash pathogens.list 0.1
kmer-db build-mh pathogens.list pathogens.db
kmer-db all2all pathogens.db matrix.csv
kmer-db distance matrix.csv
```

Calculating similarities/distances between samples listed in `pathogens.list` and `salmonella` using all k-mers. 
```
kmer-db build pathogens.list pathogens.db
kmer-db one2all pathogens.db salmonella vector.csv
kmer-db distance vector.csv
```

Same as above, but using only 10% of k-mers. 
```
kmer-db minhash pathogens.list 0.1
kmer-db build-mh pathogens.list pathogens.db
kmer-db one2all pathogens.db salmonella vector.csv
kmer-db distance vector.csv
```


