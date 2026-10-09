#include "dna_kmer.hpp"

namespace DNA
{
	KmerSet::~KmerSet()
	{
		delete[] _occurrences;
	}

	std::ostream& operator<<(std::ostream& os, KmerSet const& kmer)
	{
		for (size_t i {}; i < kmer.size(); ++i)
		{
		}

		return os;
	}
}

