#ifndef DNA_KMER_HPP
#define DNA_KMER_HPP

#include <iostream>
#include "dna_sequence.hpp"

namespace DNA 
{
	class KmerSet
	{
		public:
			constexpr KmerSet(Sequence const& dna_sequence, const size_t k)
				: _size(dna_sequence.get_kmer_count(k))
				, _occurrences(new short[1ULL<<(2*k)]())
			{
				for (size_t i {}; i < _size; ++i)
					++_occurrences[dna_sequence.get_kmer(i, k)];
			}

			~KmerSet();

			inline constexpr size_t size() const { return _size; }
			inline constexpr uint32_t operator[](const size_t i) const 
			{
				return _occurrences[i];
			}

		private:
			const size_t _size;
			short *const _occurrences;
			
			// Size: 16B 
			// Alignment: 8B
	};
}

#endif
