#include "dna_similarity.hpp"
#include "dna_kmer.hpp"
#include "dna_sequence.hpp"

/* Contant expressions */

static constexpr uint32_t max(const uint32_t a, const uint32_t b) 
{ return (a+b+std::abs((int)a-(int)b))/2; }

static constexpr uint32_t min(const uint32_t a, const uint32_t b)
{ return (a+b-std::abs((int)a-(int)b))/2; }

/* Private declarations */

/*
 * Calculates the similarities between two DNA sequences.
 * Uses a contigious array with occurrances indexed by the 
 * kmer value. 
 */
static double calculate_using_kmer_array(
		DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, 
		const size_t k);

// TODO: Add comment
static double calculate_using_kmer_buckets(
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
			return calculate_using_kmer_array(seq_a, seq_b, k);
		
		return calculate_using_kmer_buckets(seq_a, seq_b, k);
	}
}

/* Private definitions */

// TODO: Add vectorization
static double calculate_using_kmer_array(DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, const size_t k)
{
	const DNA::KmerArray kmer_arr { seq_a, seq_b, k };
	const size_t max_kmer {1ULL << (2*k)};
	const size_t mid = seq_a.size() + seq_b.size() - 2*k + 2;

	uint64_t diff {0};
	for (size_t i {}; i < max_kmer; ++i)
		diff += std::abs(kmer_arr[i]);

	return (double)(mid-diff) / (double)(mid+diff);
}

static double calculate_using_kmer_buckets(DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, const size_t k)
{
	DNA::KmerBuckets kmer_buckets { seq_a, seq_b, k };
	const size_t mid = seq_a.size() + seq_b.size() - 2*k + 2;
	
	uint64_t diff {0};
	for (size_t bucket_i {}; bucket_i < 256; ++bucket_i)
	{
		const size_t bucket_size = kmer_buckets.get_bucket_size(bucket_i);
		size_t i {};
		while (i < bucket_size)
		{
			uint64_t kmer = kmer_buckets.at(bucket_i, i).kmer_segment;
			int count_a {0};
			int count_b {0};

			while((i < bucket_size) && 
					(kmer_buckets.at(bucket_i,i).kmer_segment == kmer))
			{
				auto const& item = kmer_buckets.at(bucket_i, i);
				count_a += (!item.from_b) * item.count; 
				count_b += item.from_b * item.count;
				
				++i;
			}

			diff += std::abs(count_a - count_b);
		}
	}

	return (double)(mid-diff) / (double)(mid+diff);
}

