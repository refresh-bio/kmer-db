#include "tests.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void Tests::comparePatterns(const AbstractKmerDb& db1, const AbstractKmerDb& db2, const std::vector<kmer_t>& kmers) {

	cout << "Comparing databases..." << endl;

	int i = 0;

	std::vector<sample_id_t> samples1;
	std::vector<sample_id_t> samples2;

	for (auto kmer : kmers) {
		db1.mapKmers2Samples(kmer, samples1);
		db2.mapKmers2Samples(kmer, samples2);
			
		bool eq = std::equal(samples1.begin(), samples1.end(), samples2.begin(), samples2.end());

		if (!eq) {
			std::cout << "k-mer: " << kmer << std::endl;
		}

	}
}


void Tests::testSerialization(const FastKmerDb& db) {
	
	cout << "Serializing..." << endl;

	std::ofstream ofs("D:/kmer-db.bin", std::ios::binary);
	db.serialize(ofs);
	ofs.close();

	cout << "Deserializing..." << endl;

	std::ifstream ifs("D:/kmer-db.bin", std::ios::binary);
	FastKmerDb copy;
	copy.deserialize(ifs);
	ifs.close();

	cout << "Comparing...";
	comparePatterns(db, copy, db.getKmers());

	cout << "Done!" << endl;

}


 void Tests::testDistanceMatrix(const AbstractKmerDb& db, const string& fname) {

	 cout << "Calculating distance matrix..." << endl;

	 std::ofstream file(fname);

	 file << endl << "Distance matrix:" << endl;
	 Array<uint32_t> matrix;
	 db.calculateSimilarityMatrix(matrix);
	 cout << "Saving distance matrix...";
	 for (int i = 0; i < matrix.size(); ++i) {
		 for (int j = 0; j < matrix.size(); ++j) {
			 file << setw(10) << matrix[i][j];
		 }
		 file << endl;
	 }

	 cout << "done!" << endl;

/*	 fileNaive << endl << "Histogram: " << endl;
	 auto stats = naive_db.getPatternsStatistics();
	 for (auto s : stats) {
		 fileNaive << setw(10) << s.second << ": ";
		 copy(s.first.begin(), s.first.end(), ostream_iterator<sample_id_t>(fileNaive, ","));
		 fileNaive << endl;
	 }
	 */
}