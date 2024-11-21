#pragma once

#define VERSION "2.2.5"
#define DATE "02.11.2024"

/* 

Version history

2.2.5 (02.11.2024)
- Improved logging.

2.2.4 (31.10.2024)
- Improved support of patterns shared by large number of samples ('bubbles').

2.2.3 (21.10.2024)
- Added `-version` switch.

2.2.2 (04.10.2024)
- Fixed bug with filtering/sampling in all2all variants and new2all always using 18 as a kmer length.

2.2.1 (03.10.2024)
- Fixed slight slow down in `build` mode introduced in the previous version. 

2.2.0 (02.10.2024)
- Added rows sampling (`-sample-rows` option) in `all2all-sp` and `all2all-parts` modes. 

2.1.0 (27.09.2024)
- Added unified filtering based on any specified measure (parameters `-above`, `-above_eq`, `-below`, and `-below_eq` replaced with `-min` and `-max` options).
- Changed interface in `distance` mode: only one measure allowed, output file has to be specified.

2.0.6 (23.09.2024)
- Speed-ups in all2all-sp and all2all-parts modes.

2.0.5 (12.09.2024)
- Updates in tests and automatic building scripts.

2.0.4 (20.08.2024)
- Above and below options working correctly in all2all-sp mode.

2.0.3 (28.06.2024)
- Fixed bug with empty sample.

2.0.2 (19.06.2024)
- Some unneccessary stuff removed from the database.

2.0.1 (18.06.2024)
- Improved parallelization scheme.


2.0.0 (31.05.2024)
- Added new modes:  all2all_sparse, all2all_parts,
- Serious time and memory optimizations,
- Support of MacOS (Apple and x86 CPUs) and ARM platforms.


1.11.1 (07.03.2023)
- Removed deadlock in the -multisample-fasta mode.

1.11.0 (27.02.2023)
- Added -below and -above thresholds for all2all and new2all modes.

1.10.0 (18.07.2022)
- Added support of sparse inputs in distance mode,
- Added support of sparse outputs in distance mode (-sparse switch) with optional filtering (-above/-below options),
- Extended help information.

1.9.4 (19.04.2022)
- fixed database construction for very small samples (#kmers < #threads)
- fixed synchronization issues in new2all mode (non-deterministic row order in output matrix).
- fixed deadlock during database construction when -multisample-fasta mode is run on more than one file.

1.9.3 (27.08.2021)
- Disabled h-mer hashatables loading in all2all mode.
- Fast CSV saving in all2all and new2all modes.

1.9.2 (16.08.2021)
- Synchronization bugfix in new2all.

1.9.1 (11.08.2021)
- Output matrices can be stored in sparse format (-sparse switch).
- Better workload balancing.

1.9.0 (09.08.2021)
- Improved parallelization scheme in new2all mode (few-fold speed improvement).
- Reduced memory footprint of -multisample-fasta mode.
- More than one input FASTA files supported in -multisample-fasta mode.

1.8.0 (19.03.2021)
- Added -extend switch which allows extending existing kmer database.
- Serialization/deserialization works much faster now.
- Fixed serious bug in -multisample-fasta mode which caused incorrect kmers counting.

1.7.6 (31.03.2020)
- Fixed bug in distance mode when sequence id contained spaces.

1.7.5 (13.02.2020)
- Some compilation warnings removed.
- Fixed crash on samples with small k-mers count or very small filter values.

1.7.4 (29.01.2020)
- Proper handling of triangle input matrices in `distance` mode (triangle outputs are produced).

1.7.3 (17.01.2020)
- Fixed rare bug in hashtable when k-mer containing only T bases was treated as an empty entry. Now an empty item is indicated by a special value instead of a special key.

1.7.2 (15.01.2020)
- Added new distance measure `-mash-query` which is a mash distance calculated w.r.t. a query length (use if the query is much shorter than database sequences).
- C++11 compatibility (compiles with G++ 4.8.5).  

1.7.1 (31.10.2019)
- Possibility to specify low threshold of k-mer minhash filter (-f-start parameter).
- When loading genome files, exact filenames are examined first. If this fails, an attempt to add predefined extensions is made.  

1.7.0 (27.09.2019)
- For performance reasons upper triangle (with diagonal) of distance matrix in all2all mode is no longer saved.
- Preparations for raw serialization of hashtables.

1.6.2: (28.05.2019) After data structure update - stable
Note: Starting from this release version numbering conforms to major.minor.patch scheme.
Added:
- Switch-phylip-out in distance mode which allows storing distance/similarity matrices in Phylip format.

Fixed several bugs from 1.51 release:
- Incorrect support of k-mer lengths < 16.
- Very long processing of long k-mers (k >= 26).
- Segmentation fault when storing minhashed k-mers on a disk (minhash mode).


1.51 (11.04.2019)
- Serious reduction of time and memory requirements of build mode caused by the changes of the data structures. 
E.g., when tested on full k-mer spectrum of 40715 pathogen genomes, time and memory footprint decreased by 1/3 (1h30 to 1h, 60 to 40GB).
- Several new parameters added.
- Lots of bugs fixed.

1.20 (12.02.2019)
Changes:
- new2all mode added,
- uniform output table format for all2all, new2all, and one2all modes,

Bugs fixed:
- proper support of samples with no k-mers of given length,
- problems with building database from minhashed k-mers.

1.12 (19.07.2018)
- Support of no-AVX2 build.

1.11 (13.07.2018)
- Proper handling of uncompressed genomes. 
- Linking against pre-installed zlib available.

1.1 (2018-06-12)
- File loading refactored.
- Small bugs fixed.

1.0 (10.02.2018)
Initial release














*/
