set EXEC="./../../src/x64/Release/kmer-db.exe"
set INPUT_DIR=./protein

cd %INPUT_DIR%


%EXEC% build -t 1 -multisample-fasta -k 24 dna_100x1000.fasta dna.db
%EXEC% build -t 1 -multisample-fasta -k 24 -preserve-strand dna_100x1000.fasta dna-preserve.db 
%EXEC% build -t 1 -multisample-fasta -k 8 -alphabet aa aa_100x1000.fasta aa.db 

%EXEC% build -t 1 -multisample-fasta -k 7 -alphabet aa aa_100x1000.fasta aa_k7.db
 
%EXEC% build -t 1 -multisample-fasta -k 8 -alphabet aa11_diamond aa_100x1000.fasta aa11_diamond.db
%EXEC% build -t 1 -multisample-fasta -k 8 -alphabet aa12_mmseqs aa_100x1000.fasta aa12_mmseqs.db 
%EXEC% build -t 1 -multisample-fasta -k 8 -alphabet aa6_dayhoff aa_100x1000.fasta aa6_dayhoff.db    

rem dense inputs
%EXEC% all2all -t 1 dna.db dna.a2a
%EXEC% all2all -t 1 dna-preserve.db dna-preserve.a2a
%EXEC% all2all -t 1 aa.db aa.a2a
%EXEC% all2all -t 1 aa_k7.db aa_k7.a2a
%EXEC% all2all -t 1 aa11_diamond.db aa11_diamond.a2a
%EXEC% all2all -t 1 aa12_mmseqs.db aa12_mmseqs.a2a
%EXEC% all2all -t 1 aa6_dayhoff.db aa6_dayhoff.a2a


cd ../