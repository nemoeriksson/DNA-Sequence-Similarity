#include "dna_kmer.hpp"

namespace DNA
{
	KmerSet::~KmerSet()
	{
		delete[] _occurrences;
	}
}

