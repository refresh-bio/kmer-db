all: kmer-db

## USER'S OPTIONS
INTERNAL_ZLIB = false
NO_AVX2 = false

## ###################
KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = src
KMER_DB_LIBS_DIR = src/kmc_api

ifeq ($(INTERNAL_ZLIB),true)
	EXTRA_LIBS_DIR = libs
else
	EXTRA_LIBS_DIR = ""
endif

CC = g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -fopenmp -pthread -mavx -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CFLAGS_AVX2	= -Wall -O3 -m64 -std=c++14 -fopenmp -pthread -mavx2 -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CLINK	= -lm -O3 -std=c++14 -lpthread -fopenmp -mavx -fabi-version=6 

OBJS := $(KMER_DB_MAIN_DIR)/kmer_db.o \
	$(KMER_DB_MAIN_DIR)/analyzer.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/instrset_detect.o \
	$(KMER_DB_MAIN_DIR)/loader.o \
	$(KMER_DB_MAIN_DIR)/loader_ex.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/pattern.o \
	$(KMER_DB_MAIN_DIR)/kmc_file_wrapper.o \
	$(KMER_DB_MAIN_DIR)/prefix_kmer_db.o \
	$(KMER_DB_MAIN_DIR)/similarity_calculator.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o 

$(KMER_DB_MAIN_DIR)/parallel_sorter.o: $(KMER_DB_MAIN_DIR)/parallel_sorter.cpp
	$(CC) -O3 -mavx -m64 -std=c++14 -pthread -fopenmp -c $< -o $@

ifeq ($(NO_AVX2),true)
## no avx2 support
AVX_OBJS := $(KMER_DB_MAIN_DIR)/row_add_avx.o 
$(KMER_DB_MAIN_DIR)/row_add_avx.o: $(KMER_DB_MAIN_DIR)/row_add_avx.cpp
	$(CC) $(CFLAGS) -DNO_AVX2 -c $< -o $@

else
# with avx2 support
AVX_OBJS := $(KMER_DB_MAIN_DIR)/row_add_avx.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx2.o 
$(KMER_DB_MAIN_DIR)/row_add_avx.o: $(KMER_DB_MAIN_DIR)/row_add_avx.cpp
	$(CC) $(CFLAGS) -c $< -o $@
$(KMER_DB_MAIN_DIR)/row_add_avx2.o: $(KMER_DB_MAIN_DIR)/row_add_avx2.cpp
	$(CC) $(CFLAGS_AVX2) -c $< -o $@

endif


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

ifeq ($(INTERNAL_ZLIB),true)
kmer-db: $(OBJS) $(AVX_OBJS)
	$(CC) $(CLINK) -o $(KMER_DB_ROOT_DIR)/$@ $(OBJS) $(AVX_OBJS) $(EXTRA_LIBS_DIR)/libz.a
else
kmer-db: $(OBJS) $(AVX_OBJS)
	$(CC) $(CLINK) -o $(KMER_DB_ROOT_DIR)/$@ $(OBJS) $(AVX_OBJS) -lz
endif	

clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm kmer-db
	
