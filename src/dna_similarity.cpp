#include "dna_similarity.hpp"
#include "dna_kmer.hpp"

namespace DNA 
{
	static constexpr uint32_t max(const uint32_t a, const uint32_t b) 
	{
		return a > b ? a : b;
	}

	static constexpr uint32_t min(const uint32_t a, const uint32_t b)
	{
		return a < b ? a : b;
	}

	double get_similarity(Sequence const& seq_a, 
			Sequence const& seq_b, 
			const size_t k)
	{
		const KmerSet set_a { seq_a, k };
		const KmerSet set_b { seq_b, k };
		const size_t max_kmer {1ULL << (2*k)};

		double intersection_sum {};
		double union_sum {};

		for (size_t i {}; i < max_kmer; ++i)
		{
			intersection_sum += min(set_a[i], set_b[i]);
			union_sum += max(set_a[i], set_b[i]);
		}

		return intersection_sum / union_sum;
	}
}

