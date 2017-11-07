all: kmer_db

KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = kmer_db
KMER_DB_LIBS_DIR = kmer_db/kmc_api

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++11 -fopenmp -pthread -static -I $(KMER_DB_LIBS_DIR)
CLINK	= -lm -O3 -std=c++11 -lpthread -fabi-version=6 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

kmer_db: $(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/queue.o \
	$(KMER_DB_MAIN_DIR)/tests.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o
	$(CC) $(CLINK) -o $(KMER_DB_ROOT_DIR)/$@  \
	$(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/queue.o \
	$(KMER_DB_MAIN_DIR)/tests.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o

clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm kmer_db
	
