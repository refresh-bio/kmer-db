all: kmer-db-1.0

KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = src
KMER_DB_LIBS_DIR = src/kmc_api
EXTRA_LIBS_DIR = libs

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -fopenmp -pthread -mavx -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CFLAGS_AVX2	= -Wall -O3 -m64 -std=c++14 -fopenmp -pthread -mavx2 -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CLINK	= -lm -O3 -std=c++14 -lpthread -fopenmp -mavx -fabi-version=6 

$(KMER_DB_MAIN_DIR)/parallel_sorter.o: $(KMER_DB_MAIN_DIR)/parallel_sorter.cpp
	$(CC) -O3 -mavx -m64 -std=c++14 -pthread -fopenmp -c $< -o $@

$(KMER_DB_MAIN_DIR)/row_add_avx2.o: $(KMER_DB_MAIN_DIR)/row_add_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

kmer-db-1.0: $(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/instrset_detect.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx2.o \
	$(KMER_DB_MAIN_DIR)/loader.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/pattern.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o 
	$(CC) $(CLINK) -o $(KMER_DB_ROOT_DIR)/$@  \
	$(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/instrset_detect.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx2.o \
	$(KMER_DB_MAIN_DIR)/loader.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/pattern.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o \
	$(EXTRA_LIBS_DIR)/libz.a \
	$(EXTRA_LIBS_DIR)/libbz2.a \
	
clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm kmer-db-1.0
	
