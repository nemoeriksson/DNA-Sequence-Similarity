#include "dna_similarity.hpp"
#include "dna_kmer.hpp"
#include "dna_sequence.hpp"

/* Contant expressions */

static constexpr uint32_t max(const uint32_t a, const uint32_t b) 
{ return a > b ? a : b; }

static constexpr uint32_t min(const uint32_t a, const uint32_t b)
{ return a < b ? a : b; }

/* Private declarations */

/*
 * Calculates the similarities between two DNA sequences.
 * Uses a contigious array with occurrances indexed by the 
 * kmer value. 
 */
static double calculate_using_kmer_set(
		DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, 
		const size_t k);

/* Public definitions */

namespace DNA 
{
	double get_similarity(Sequence const& seq_a, 
			Sequence const& seq_b, 
			const size_t k)
	{
		if (k <= 8)
			return calculate_using_kmer_set(seq_a, seq_b, k);
	
		return -1;
	}
}

/* Private definitions */

// TODO: Add vectorization
static double calculate_using_kmer_set(DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, const size_t k)
{
	const DNA::KmerSet set_a { seq_a, k };
	const DNA::KmerSet set_b { seq_b, k };
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

