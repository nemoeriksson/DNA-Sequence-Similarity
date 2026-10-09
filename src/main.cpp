#include <iostream>
#include "dna_sequence.hpp"
#include "dna_kmer.hpp"

int main(void)
{
	DNA::Sequence dna_seq {"ACTG"};
	std::cout << dna_seq << std::endl;
	DNA::KmerSet kmer_set {dna_seq, 2};

	return 0;
}

