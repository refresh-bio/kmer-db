all: kmer-db-1.0

KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = kmer_db
KMER_DB_LIBS_DIR = kmer_db/kmc_api
KMER_DB_RADULS_DIR = kmer_db/raduls

CC 	= /usr/local/gcc62/bin/g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -fopenmp -pthread -I $(KMER_DB_LIBS_DIR)
CLINK	= -lm -O3 -std=c++14 -lpthread -fopenmp -fabi-version=6 

$(KMER_DB_RADULS_DIR)/sorting_network.o: $(KMER_DB_RADULS_DIR)/sorting_network.cpp
	$(CC) -O1 -m64 -std=c++14 -mavx -pthread -c $< -o $@

$(KMER_DB_MAIN_DIR)/parallel_sorter.o: $(KMER_DB_MAIN_DIR)/parallel_sorter.cpp
	$(CC) -O3 -mavx -m64 -std=c++14 -fopenmp -pthread -fno-ipa-ra -fno-tree-vrp -fno-tree-pre -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

kmer-db-1.0: $(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/loader.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/queue.o \
	$(KMER_DB_MAIN_DIR)/tests.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o \
	$(KMER_DB_RADULS_DIR)/sorting_network.o
	$(CC) $(CLINK) -o $(KMER_DB_ROOT_DIR)/$@  \
	$(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/loader.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/queue.o \
	$(KMER_DB_MAIN_DIR)/tests.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o \
	$(KMER_DB_RADULS_DIR)/sorting_network.o
	

clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm $(KMER_DB_RADULS_DIR)/*.o
	-rm kmer-db-1.0
	
