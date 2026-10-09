#ifndef DNA_KMER
#define DNA_KMER

#include <iostream>
#include "dna_sequence.hpp"

namespace DNA 
{
	class KmerSet
	{
		public:
			constexpr KmerSet(Sequence const& dna_sequence, const size_t k)
				: _size(dna_sequence.get_kmer_count(k))
				, _occurrences(new short[_size]())
			{
				for (size_t i {}; i < _size; ++i)
				{
					uint32_t kmer {};
					for (size_t j {}; j < k; ++j)
					{
						kmer <<= 2;
						kmer |= dna_sequence[i+j];
						std::cout << Sequence::to_letter(dna_sequence[j+i]);
					}
					std::cout << std::endl;

					++_occurrences[kmer];
				}
			}

			~KmerSet();

			constexpr inline size_t size() const { return _size; }


		private:
			const size_t _size;
			short *const _occurrences;
			
			// Size: 16B 
			// Alignment: 8B
	};
}

#endif
