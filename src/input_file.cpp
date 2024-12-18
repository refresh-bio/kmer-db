/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include "input_file.h"
#include "kmc_api/kmc_file.h"
#include "kmer_db.h"
#include "filter.h"
#include "kmer_extract.h"
#include "parallel_sorter.h"

#include "../libs/refresh/compression/lib/file_wrapper.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>
#include <filesystem>

#ifdef USE_RADULS
	#include <raduls.h>
#endif

using namespace std;

