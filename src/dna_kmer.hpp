#ifndef DNA_KMER_HPP
#define DNA_KMER_HPP

#include <iostream>
#include "dna_sequence.hpp"

namespace DNA 
{
	/*
	 * Stores the occurrences of each k-mer value in a 
	 * contigious array, indexed by the k-mer values.
	 *
	 * ===== Stats & Storage Usage (*_occurrences) =====
	 * Size = 2 * 4^k bytes
	 *
	 * My device:
	 * Fills 100% of L1 cache at k = 7	| 2 * 4^7  =  32 KiB of  32 KiB
	 *		  50%    L2			 k = 9	| 2 * 4^9  = 0.5 MiB of   1 MiB
	 *		  50%	 L3			 k = 11 | 2 * 4^11 =   8 MiB of  16 MiB
	 *
	 * University device (salt):
	 * Fills 100% of L1 cache at k = 7	| 2 * 4^7  =  32 KiB of  32 KiB
	 *		  50%    L2			 k = 8	| 2 * 4^8  = 128 KiB of 256 KiB
	 *		  50%	 L3			 k = 11 | 2 * 4^11 =   8 MiB of  16 MiB
	 *
	 * Recommended breakpoint: k=8 
	 */
	class KmerSet
	{
		public:
			KmerSet(Sequence const& dna_sequence, const size_t k);
			
			~KmerSet();

			inline size_t size() const { return _size; }

			inline uint16_t operator[](const size_t i) const 
			{
				return _occurrences[i];
			}

		private:
			const size_t _size;
			uint16_t *const _occurrences;
			
			// Size: 16B 
			// Alignment: 8B
	};
}

#endif
