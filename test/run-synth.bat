set EXEC="./../../src/x64/Release/kmer-db.exe"
set INPUT_DIR=./synth

cd %INPUT_DIR%

#%EXEC% build -multisample-fasta -k 21 synth-local.list synth.db

rem dense inputs
%EXEC% all2all synth.db a2a
#%EXEC% distance mash ani a2a
#%EXEC% distance -sparse ani a2a
#%EXEC% distance -sparse -above -1.0 -below 1.0 mash a2a

#%EXEC% new2all -multisample-fasta synth.db synth-local.list n2a
#%EXEC% distance mash ani n2a
#%EXEC% distance -sparse mash ani n2a

rem sparse inputs
#%EXEC% all2all -sparse synth.db a2a-sparse
#%EXEC% distance mash ani a2a-sparse
#%EXEC% distance -sparse mash ani a2a-sparse

#%EXEC% new2all -multisample-fasta -sparse synth.db synth-local.list n2a-sparse
#%EXEC% distance mash ani n2a-sparse
#%EXEC% distance -sparse mash ani n2a.sparse.csv

%EXEC% all2all-sp synth.db a2a-sp

cd ../