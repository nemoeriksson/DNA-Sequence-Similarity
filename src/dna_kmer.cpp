#include "dna_kmer.hpp"

namespace DNA
{			
	KmerSet::KmerSet::KmerSet(Sequence const& dna_sequence, const size_t k)
		: _size(dna_sequence.get_kmer_count(k))
		, _occurrences(new short[1ULL<<(2*k)]())
	{
		for (size_t i {}; i < _size; ++i)
			++_occurrences[dna_sequence.get_kmer(i, k)];
	}


	KmerSet::~KmerSet()
	{
		delete[] _occurrences;
	}
}

