set EXEC="./../../src/x64/Release/kmer-db.exe"
set INPUT_DIR=./virus

rem %EXEC% build %INPUT_DIR%/seqs.list k18.db

rem %EXEC% all2all k18.db k18.csv

rem %EXEC% all2all -sparse k18.db k18.sparse.csv

rem %EXEC% build -f 0.1 %INPUT_DIR%/seqs.list k18.frac.db
rem %EXEC% all2all k18.frac.db k18.frac.csv

rem %EXEC% minhash 0.1 %INPUT_DIR%/seqs.list
rem %EXEC% build -from-minhash -k 25 %INPUT_DIR%/seqs.list k18.minhash.db
rem %EXEC% all2all k18.minhash.db k18.minhash.csv


rem %EXEC% build -multisample-fasta %INPUT_DIR%/multi.list k18.multi.db
rem %EXEC% all2all k18.multi.db k18.multi.csv

rem %EXEC% build %INPUT_DIR%/seqs.part1.list k18.parts.db
rem %EXEC% build -extend -k 25 %INPUT_DIR%/seqs.part2.list k18.parts.db
rem %EXEC% all2all k18.parts.db k18.parts.csv
		
rem %EXEC% build -k 24 %INPUT_DIR%/seqs.list k24.db
rem %EXEC% all2all k24.db k24.csv

rem %EXEC% build %INPUT_DIR%/seqs.part1.list k18.parts.22.db
rem %EXEC% new2all k18.parts.22.db %INPUT_DIR%/seqs.part2.list k18.n2a.22.csv

rem %EXEC% new2all -sparse k18.parts.db %INPUT_DIR%/seqs.part2.list k18.n2a.sparse.csv

rem %EXEC% new2all k18.db %INPUT_DIR%/seqs.list k18.n2a.itself.csv
		
rem %EXEC% distance jaccard min max cosine mash k18.csv


rem %EXEC% build -k 25 -f 0.1 -t 16  %INPUT_DIR%/seqs.part1.list k25.db
rem %EXEC% one2all k25.db %INPUT_DIR%/data/MT159713 MT159713.csv

cd %INPUT_DIR%

%EXEC% build seqs.part1-local.list k18.part1.db
%EXEC% build seqs.part2-local.list k18.part2.db
%EXEC% all2all-parts multi.db.list k18.parts.csv

cd ..