/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#ifdef _MSC_VER 
//#include <mimalloc.h>
#include <mimalloc-new-delete.h>
#endif

#include "log.h"
#include "version.h"

#include "console.h"
#include "params.h"

#include <stdexcept>

int main(int argc, char **argv)
{
	Log::getInstance(Log::LEVEL_NORMAL).enable();
	
	Params params;

	try {
		// returns false when help message was desired
		if (!params.parse(argc, argv)) {
			return 0;
		}

		time_t rawtime;
		struct tm* timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		LOG_NORMAL("Analysis started at " << asctime(timeinfo) << endl);

		auto console = ConsoleFactory::create(params.mode);
		if (!console) {
			throw std::runtime_error("Invalid mode selected");
		}

		console->run(params);

		time(&rawtime);
		timeinfo = localtime(&rawtime);
		LOG_NORMAL(endl << "Analysis finished at " << asctime(timeinfo) << endl);
	}
	catch (usage_error& err) {
		LOG_NORMAL("ERROR: Incorrect usage" << endl << "See detailed instructions below" << endl << endl);
		params.showInstructions(err.getMode());
		return -1;
	}
	catch (std::runtime_error& err) {
		LOG_NORMAL("ERROR: " << err.what() << endl);
		return -1;
	}

	return 0;
}
