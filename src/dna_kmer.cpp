#include "dna_kmer.hpp"

namespace DNA
{			
	KmerSet::KmerSet::KmerSet(Sequence const& dna_sequence, const size_t k)
		: _size(dna_sequence.get_kmer_count(k))
		, _occurrences(new uint16_t[1ULL<<(2*k)]())
	{
		uint64_t kmer {};
		for (size_t i {}; i < k; ++i)
			kmer = (kmer<<2)|dna_sequence[i];
		for (size_t i {}; i < _size; ++i)
		{
			kmer = ((kmer<<2)|dna_sequence[i]) & ((1<<(2*k))-1);
			++_occurrences[kmer];
		}
	}


	KmerSet::~KmerSet()
	{
		delete[] _occurrences;
	}
}

