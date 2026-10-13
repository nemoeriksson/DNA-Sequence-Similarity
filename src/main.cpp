#include <chrono>
#include <cstdlib>
#include <iostream>
#include "dna_sequence.hpp"
#include "dna_similarity.hpp"

int main()
{
	constexpr size_t iterations = 1000;
	constexpr size_t max_k = 8;

	DNA::Sequence seq_a {"ATACGATCGATATCTAGCTAGCTAGCTATATCTAGCATTAGCTAGCTATATCGATCGATATCGATGCTATAGCTAGTATTATATATAAGATGCATATTATAAAAATATATTATATCGATATTACGATGCTAGCTATATATCGATCGTATATCATTATATCGATATAT"};
	DNA::Sequence seq_b {"AAATCATCATGTCAGCTACTATCTACTGACTACTACTAGCATTACATAATTCATGATCATCGAGCTACAGATTATATCGATATCTAGCTAGCTATATTATCATGCTAGCTATATATTATATCGATCAGTATATATAATATAGCGCTAATAATATCGATATCGATGCT"};
	
	double similarity;

	std::cout << "k;Similarity Score; Average Calculation Time" << std::endl; // CSV header

	for (size_t k {1}; k <= max_k; ++k)
	{
		auto start_time = std::chrono::system_clock::now();
		for (size_t i {}; i < iterations; ++i)
			similarity = DNA::get_similarity(seq_a, seq_b, k);
		auto end_time = std::chrono::system_clock::now();

		auto elapsed = end_time - start_time;

		// CSV data
		std::cout << k << ';';
		std::cout << similarity << ';';
		std::cout << elapsed/iterations;
		std::cout << std::endl;
	}

	return 0;
}

