#ifndef DNA_KMER_HPP
#define DNA_KMER_HPP

#include <iostream>
#include <vector>
#include "dna_sequence.hpp"

namespace DNA 
{
	class KmerArray
	{
		using difference_type = int16_t;

		public:
			KmerArray(Sequence const& seq_a, 
					Sequence const& seq_b,
					const size_t k);
			
			~KmerArray();

			inline difference_type operator[](const size_t i) const 
			{
				return _differences[i];
			}

		private:
			difference_type *const _differences;
	};

	class KmerBuckets
	{
		public:
			struct BucketItem {
				uint64_t from_b: 1;
				uint64_t count: 7;
				uint64_t kmer_segment: 56;
			}; 

			KmerBuckets(Sequence const& seq_a,
				Sequence const& seq_b,
				const size_t k);

			inline size_t get_bucket_size(const size_t bucket_index) const
			{
				return _buckets[bucket_index].size();
			}

			BucketItem const& at(const size_t bucket_index, const size_t i) const 
			{
				return _buckets[bucket_index][i];
			}

		private:
			std::vector<BucketItem> _buckets[256] {};
	};
}

#endif
