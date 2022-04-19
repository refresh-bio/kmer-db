all: kmer-db

## USER'S OPTIONS
INTERNAL_ZLIB = false

####################

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

ifeq ($(uname_S),Linux)   
	# check if CPU supports AVX2
	HAVE_AVX2=$(filter-out 0,$(shell grep avx2 /proc/cpuinfo | wc -l))
	OMP_FLAGS = -fopenmp
	ABI_FLAGS = -fabi-version=6
endif
ifeq ($(uname_S),Darwin)
	 # check if CPU supports SSE4.2
	HAVE_AVX2=$(filter-out 0,$(shell  sysctl -n machdep.cpu.features machdep.cpu.leaf7_features| grep AVX2 - | wc -l))
	OMP_FLAGS = -Xpreprocessor -fopenmp 
	ABI_FLAGS = 
endif

## ###################
KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = src
KMER_DB_LIBS_DIR = src/kmc_api

ifeq ($(INTERNAL_ZLIB),true)
	EXTRA_LIBS_DIR = libs
else
	EXTRA_LIBS_DIR = ""
endif

LDFLAGS += 
CFLAGS	+= -Wall -O3 -m64 -std=c++11 $(OMP_FLAGS) -pthread 
CFLAGS_AVX2	+= $(CFLAGS) -mavx2 -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CFLAGS += -mavx  -I $(KMER_DB_LIBS_DIR) -I $(EXTRA_LIBS_DIR)
CLINK	= -lm -O3 -std=c++11 -lpthread $(OMP_FLAGS) -mavx $(ABI_FLAGS) 


OBJS := $(KMER_DB_MAIN_DIR)/analyzer.o \
	$(KMER_DB_MAIN_DIR)/console.o \
	$(KMER_DB_MAIN_DIR)/instrset_detect.o \
	$(KMER_DB_MAIN_DIR)/loader_ex.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/pattern.o \
	$(KMER_DB_MAIN_DIR)/input_file.o \
	$(KMER_DB_MAIN_DIR)/prefix_kmer_db.o \
	$(KMER_DB_MAIN_DIR)/similarity_calculator.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o 

$(KMER_DB_MAIN_DIR)/parallel_sorter.o: $(KMER_DB_MAIN_DIR)/parallel_sorter.cpp
	$(CXX) -O3 -mavx -m64 -std=c++11 -pthread $(OMP_FLAGS) -c $< -o $@

ifeq ($(HAVE_AVX2),)
## no avx2 support
AVX_OBJS := $(KMER_DB_MAIN_DIR)/row_add_avx.o 
$(KMER_DB_MAIN_DIR)/row_add_avx.o: $(KMER_DB_MAIN_DIR)/row_add_avx.cpp
	$(CXX) $(CFLAGS) -DNO_AVX2 -c $< -o $@

else
# with avx2 support
AVX_OBJS := $(KMER_DB_MAIN_DIR)/row_add_avx.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx2.o 
$(KMER_DB_MAIN_DIR)/row_add_avx.o: $(KMER_DB_MAIN_DIR)/row_add_avx.cpp
	$(CXX) $(CFLAGS) -c $< -o $@
$(KMER_DB_MAIN_DIR)/row_add_avx2.o: $(KMER_DB_MAIN_DIR)/row_add_avx2.cpp
	$(CXX) $(CFLAGS_AVX2) -c $< -o $@

endif


%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

ifeq ($(INTERNAL_ZLIB),true)
kmer-db: $(OBJS) $(AVX_OBJS)
	$(CXX) $(CLINK) $(LDFLAGS) -o $(KMER_DB_ROOT_DIR)/$@ $(OBJS) $(AVX_OBJS) $(EXTRA_LIBS_DIR)/libz.a
else
kmer-db: $(OBJS) $(AVX_OBJS)
	$(CXX) $(CLINK) $(LDFLAGS) -o $(KMER_DB_ROOT_DIR)/$@ $(OBJS) $(AVX_OBJS) -lz 
endif	

clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm kmer-db
	
