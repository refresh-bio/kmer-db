all: kmer-db

ifeq ($(FORCE_ZLIB),true)
$(info Forcing ZLIB-NG linking)
else
NASM_V := $(shell nasm --version 2>/dev/null)
endif

####################

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
    uname_M := "x86_64"
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
    uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

ifeq ($(uname_S),Linux)   
	# check if CPU supports AVX2
	HAVE_AVX2=$(filter-out 0,$(shell grep avx2 /proc/cpuinfo | wc -l))
	ABI_FLAGS = -fabi-version=6
endif
ifeq ($(uname_S),Darwin)
	 # check if CPU supports SSE4.2
	HAVE_AVX2=$(filter-out 0,$(shell  sysctl -n machdep.cpu.features machdep.cpu.leaf7_features| grep AVX2 - | wc -l))
	ABI_FLAGS = 
endif

ifeq ($(PLATFORM), arm8)
$(info *** ARMv8 with NEON extensions ***)
	ARCH_FLAGS := -march=armv8-a  -DARCH_ARM
else ifeq ($(PLATFORM), m1)
$(info *** Apple M1(or never) with NEON extensions ***)
	ARCH_FLAGS := -march=armv8.4-a  -DARCH_ARM
else ifeq ($(PLATFORM), sse2)
$(info *** x86-64 with SSE2 extensions ***)
	ARCH_FLAGS := -msse2 -m64 -DARCH_X64 
else ifeq ($(PLATFORM), avx)
$(info *** x86-64 with AVX extensions ***)
	ARCH_FLAGS := -mavx -m64  -DARCH_X64
else ifeq ($(PLATFORM), avx2)
$(info *** x86-64 with AVX2 extensions ***)
	ARCH_FLAGS := -mavx2 -m64  -DARCH_X64
else
$(info *** Unspecified platform - use native compilation)
	ifeq ($(uname_M),x86_64)
		ARCH_FLAGS := -march=native -DARCH_X64
	else
		ARCH_FLAGS := -march=native -DARCH_ARM
	endif	
endif

CFLAGS = $(ARCH_FLAGS)

ifeq ($(uname_M),x86_64)
	ifdef NASM_V
$(info NASM detected - using Intel ISA-L ZLIB implementation)
		GZ_LIB:=isa-l.a
		GZ_INCLUDE:=libs/isa-l/include
		gz_target:=isa-l
		CFLAGS+=-DREFRESH_USE_IGZIP
	else
		GZ_LIB:=libz.a
		GZ_INCLUDE:=libs/zlib-ng/include
		gz_target:=ng_zlib
		CFLAGS+=-DREFRESH_USE_ZLIB
	endif
else
	GZ_LIB:=libz.a
	gz_target:=ng_zlib
	CFLAGS+=-DREFRESH_USE_ZLIB
endif

## ###################
KMER_DB_ROOT_DIR = .
KMER_DB_MAIN_DIR = src
KMER_DB_LIBS_DIR = src/kmc_api

MIMALLOC_OBJ=libs/mimalloc/mimalloc.o

LIBS_DIR = libs
ISAL_DIR = libs/isa-l
ZLIB_DIR = libs/zlib-ng

INC_DIRS =. libs/mimalloc/include $(GZ_INCLUDE) libs
INCLUDE_DIR=$(foreach d, $(INC_DIRS), -I$d)

CFLAGS	+= -Wall -O3 -std=c++20 -static -pthread -g
CFLAGS_AVX2	+= $(CFLAGS) -mavx2 -I $(KMER_DB_LIBS_DIR) $(INCLUDE_DIR)
CFLAGS += -I $(KMER_DB_LIBS_DIR) $(INCLUDE_DIR)
CLINK	= -lm -O3 -std=c++20 $(ABI_FLAGS) 

ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	CLINK+=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif

ifeq ($(uname_S),Darwin)
	CLINK+= -lc -static-libgcc
endif



OBJS := \
	$(KMER_DB_MAIN_DIR)/console_all2all.o \
	$(KMER_DB_MAIN_DIR)/console_all2all_parts.o \
	$(KMER_DB_MAIN_DIR)/console_all2all_sparse.o \
	$(KMER_DB_MAIN_DIR)/console_build.o \
	$(KMER_DB_MAIN_DIR)/console_db2db.o \
	$(KMER_DB_MAIN_DIR)/console_distance.o \
	$(KMER_DB_MAIN_DIR)/console_minhash.o \
	$(KMER_DB_MAIN_DIR)/console_new2all.o \
	$(KMER_DB_MAIN_DIR)/console_one2all.o \
	$(KMER_DB_MAIN_DIR)/loader_ex.o \
	$(KMER_DB_MAIN_DIR)/log.o \
	$(KMER_DB_MAIN_DIR)/main.o \
	$(KMER_DB_MAIN_DIR)/parallel_sorter.o \
	$(KMER_DB_MAIN_DIR)/params.o \
	$(KMER_DB_MAIN_DIR)/pattern.o \
	$(KMER_DB_MAIN_DIR)/input_file.o \
	$(KMER_DB_MAIN_DIR)/prefix_kmer_db.o \
	$(KMER_DB_MAIN_DIR)/similarity_calculator.o \
	$(KMER_DB_LIBS_DIR)/kmc_file.o \
	$(KMER_DB_LIBS_DIR)/kmer_api.o \
	$(KMER_DB_LIBS_DIR)/mmer.o 

$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -ftls-model=initial-exec -fno-builtin-malloc -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)


ng_zlib:
	cd $(ZLIB_DIR) && ./configure --zlib-compat && $(MAKE)  libz.a
	cp $(ZLIB_DIR)/libz.* $(LIBS_DIR)

isa-l:
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx
	cp $(ISAL_DIR)/bin/isa-l.a $(LIBS_DIR)
	cp $(ISAL_DIR)/bin/libisal.* $(LIBS_DIR)

$(KMER_DB_MAIN_DIR)/parallel_sorter.o: $(KMER_DB_MAIN_DIR)/parallel_sorter.cpp
	$(CXX) -O3 $(ARCH_FLAGS) -std=c++20 $(INCLUDE_DIR) -pthread -c $< -o $@

ifeq ($(uname_M),x86_64)
SIMD_OBJS := $(KMER_DB_MAIN_DIR)/row_add_avx.o \
	$(KMER_DB_MAIN_DIR)/row_add_avx2.o 
$(KMER_DB_MAIN_DIR)/row_add_avx.o: $(KMER_DB_MAIN_DIR)/row_add_avx.cpp
	$(CXX) $(CFLAGS) -c $< -o $@
$(KMER_DB_MAIN_DIR)/row_add_avx2.o: $(KMER_DB_MAIN_DIR)/row_add_avx2.cpp
	$(CXX) $(CFLAGS_AVX2) -c $< -o $@
else
SIMD_OBJS := $(KMER_DB_MAIN_DIR)/row_add_neon.o
$(KMER_DB_MAIN_DIR)/row_add_neon.o: $(KMER_DB_MAIN_DIR)/row_add_neon.cpp
	$(CXX) $(CFLAGS) -c $< -o $@
endif


%.o: %.cpp $(gz_target)
	$(CXX) $(CFLAGS) -c $< -o $@


kmer-db: $(OBJS) $(SIMD_OBJS) $(MIMALLOC_OBJ)
	$(CXX) -o $(KMER_DB_ROOT_DIR)/$@ $(MIMALLOC_OBJ) $(OBJS) $(SIMD_OBJS) $(LIBS_DIR)/$(GZ_LIB) $(CLINK) 


clean:
	-rm $(KMER_DB_MAIN_DIR)/*.o
	-rm $(KMER_DB_LIBS_DIR)/*.o
	-rm kmer-db
	-rm $(LIBS_DIR)/libz.*
	-rm $(LIBS_DIR)/isa-l.*
	-rm $(LIBS_DIR)/libisal.*
	-rm $(MIMALLOC_OBJ)
	cd $(ZLIB_DIR) && $(MAKE) -f Makefile.in clean
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx clean
	
