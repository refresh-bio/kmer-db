# kmer-db-dev

## USAGE

### Min-hashing k-mers
This is an optional analysis step - it can be skipped if one wants to use all k-mers from samples for distance estimation. Syntax:

`kmer-db minhash <sample_list> <fraction>`

Parameters:
 * `sample_list` - input file containing list of [KMC-generated](https://github.com/refresh-bio/KMC) k-mer files for samples to be analyzed. In the list file there should be a single entry for each sample:
    ```
    sample1
    sample2
    sample3
    ```
    When reading k-mers, the tool will automatically add `.pre` and `.suf` extensions for each sample. 
 * `fraction` - a number from <0,1> interval determining a fraction of all k-mers to be accepted by filter.
 
  For each sample from the list, a binary file with `.minhash` extension containing filtered k-mers is created.

### Building the database from k-mers

`kmer-db build <sample_list> <database>`

Parameters:
 * `sample_list` - input file containing list of k-mer files for each sample. Both, raw [KMC-generated](https://github.com/refresh-bio/KMC) as well as min-hashed k-mers are accepted. When both types of k-mer files exist, `.minhash` have higher priority.
 * `database` - output file with k-mer database for all samples.
 
 ### Calculating number of common k-mers ###
 
 `kmer_db all2all <database> <common_matrix>`
 
Parameters:
* `database`- k-mer database file

