set EXEC="./../../src/x64/Release/kmer-db.exe"
set INPUT_DIR=./ictv

cd %INPUT_DIR%

#%EXEC% build -multisample-fasta -k 25 ictv.list k25.db
#%EXEC% all2all-sp -above 10 k25.db a2a-above10.csv
#%EXEC% distance ani-shorter -above 0.7 a2a-above10.csv

#%EXEC% build -multisample-fasta -k 25 ictv.list k25.db
%EXEC% all2all-sp -min num-kmers:11 k25.db a2a-min11.csv
%EXEC% distance ani-shorter -min ani-shorter:0.7 a2a-min11.csv a2a-min11.dist-min0p7.csv

%EXEC% all2all-sp -min num-kmers:11 -min ani-shorter:0.7 k25.db a2a-min11-min0p7.csv
%EXEC% distance ani-shorter a2a-min11-min0p7.csv a2a-min11-min0p7.dist.csv
%EXEC% distance ani-shorter -min ani-shorter:0.7 a2a-min11-min0p7.csv a2a-min11-min0p7.dist-min0p7.csv

cd ..