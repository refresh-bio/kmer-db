set EXEC="./../../src/x64/Release/kmer-db.exe"
set INPUT_DIR=./synth

cd %INPUT_DIR%

%EXEC% build -multisample-fasta -k 21 synth-local.list synth.db

rem dense inputs
%EXEC% all2all synth.db a2a
%EXEC% distance mash a2a a2a.mash
%EXEC% distance ani a2a a2a.ani
%EXEC% distance -sparse ani a2a a2a.sparse-ani
%EXEC% distance -sparse -above -1.0 -below 1.0 mash a2a a2a.sparse-mash
%EXEC% distance -sparse -min 0.03 -max mash:1.0 mash a2a a2a.sparse-mash-minmax
%EXEC% distance -sparse -min 0.03 -max mash:1.0 -min num-kmers:36 mash a2a a2a.mash-sparse-min2max

%EXEC% new2all -multisample-fasta synth.db synth-local.list n2a
%EXEC% distance mash n2a n2a.mash
%EXEC% distance ani n2a n2a.ani
%EXEC% distance -sparse mash n2a n2a.sparse-mash
%EXEC% distance -sparse ani n2a n2a.sparse-ani

rem sparse inputs
%EXEC% all2all -sparse synth.db a2a-sparse
%EXEC% distance mash a2a-sparse a2a-sparse.mash
%EXEC% distance ani a2a-sparse a2a-sparse.ani
%EXEC% distance -sparse ani a2a-sparse a2a-sparse.sparse-ani
%EXEC% distance -sparse mash a2a-sparse a2a-sparse.sparse-mash

%EXEC% new2all -multisample-fasta -sparse synth.db synth-local.list n2a-sparse
%EXEC% distance mash n2a-sparse n2a-sparse.mash
%EXEC% distance ani n2a-sparse n2a-sparse.ani
%EXEC% distance -sparse mash n2a-sparse n2a-sparse.sparse-mash
%EXEC% distance -sparse ani n2a-sparse n2a-sparse.sparse-ani

%EXEC% all2all-sp synth.db a2a-sp

cd ../