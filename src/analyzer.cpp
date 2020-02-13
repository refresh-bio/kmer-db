#include "analyzer.h"

#include <iomanip>


std::unique_ptr<SetFilter> Analyzer::selectSeedKmers(const PrefixKmerDb& db, size_t minSamplesCount) {

/*	double seedFraction = (double)minSamplesCount / db.getSamplesCount();

	auto seedSet = make_unique<SetFilter>(db.getKmers2patternIds().get_size() * seedFraction); 
	const auto& kmers2patternIds = db.getKmers2patternIds();

	for (auto it = kmers2patternIds.cbegin(); it < kmers2patternIds.cend(); ++it) {
		if (kmers2patternIds.is_free(*it)) {
			continue;
		}
		pattern_id_t pid = it->val;
		const pattern_t& p = db.getPatterns()[pid];
		if (p.get_num_samples() > minSamplesCount) {
			seedSet->add(it->key);
		}
	}


	return seedSet;
	*/

	return nullptr;
}

void Analyzer::operator()(const PrefixKmerDb & db) {

}

void Analyzer::printStats(const PrefixKmerDb & db)
{

	std::vector<size_t> histo = { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, 35000, 38000, 39000  };
	std::vector<size_t> counts(histo.size());
	std::vector<size_t> kmersCounts(histo.size());

	for (const auto& p : db.getPatterns()) {
		for (size_t ih = 0; ih < histo.size(); ++ih) {
			if (p.get_num_samples() <= histo[ih]) {
				counts[ih]++;
				kmersCounts[ih] += p.get_num_kmers();
				break;
			}
		}
	}

	string desc1 = "samples count";
	string desc2 = "patterns count";
	string desc3 = "kmers count";

	cout << desc1 << ": " << desc2 << ", " << desc3 << endl;
	for (size_t ih = 0; ih < histo.size(); ++ih) {
		cout << setw(desc1.length()) << histo[ih] 
			<< ": " << Log::formatLargeNumber(counts[ih], desc2.length()) 
			<< ", " << Log::formatLargeNumber(kmersCounts[ih], desc3.length()) << endl;
	}

	//size_t minSamplesCount = 10;

	for (const auto& ht : db.getHashtables()) {
		std::vector<size_t> kmerInstances(histo.size());
	//	std::vector<size_t> uniqueKmerInstances(histo.size());

		for (auto it = ht.cbegin(); it < ht.cend(); ++it) {
			if (ht.is_free(*it)) {
				continue;
			}
			pattern_id_t pid = it->val;
			const pattern_t& p = db.getPatterns()[pid];

			for (size_t ih = 0; ih < histo.size(); ++ih) {
				if (p.get_num_samples() <= histo[ih]) {
					kmerInstances[ih] += p.get_num_samples();
					//				if (db.getRepeatedKmers().cfind(it->key) == nullptr) {
					//					uniqueKmerInstances[ih] += p.get_num_samples();
					//				}
					break;
				}
			}
		}
	}

	/*
	cout << endl << desc1 << " k-mer instances (uniques)" << endl;
	for (int ih = 0; ih < histo.size(); ++ih) {
		cout << setw(desc1.length()) << histo[ih]
			<< ": " << Log::formatLargeNumber(kmerInstances[ih])
			<< "\t(" << Log::formatLargeNumber(uniqueKmerInstances[ih]) << ")" << endl;
	}

	*/
}

