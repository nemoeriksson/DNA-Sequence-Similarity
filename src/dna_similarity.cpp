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
 * Uses an array with one entry per possible k-mer (length = 4^k)
 * where the n:th entry is equal to abs(C_A(x) - C_B(x)) where 
 * C_A(x) The count of occurrences of the k-mer x inside DNA sequence A
 * C_B(x) The count of occurrences of the k-mer x inside DNA sequence B
 *
 * Recommended for lower values k since requierd memory quickly increases
 */
static double calculate_using_kmer_array(
		DNA::Sequence const& seq_a,
		DNA::Sequence const& seq_b, 
		const size_t k);

/*
 * Splits k-mers into different "buckets" depending on their 8 most 
 * significant bits. Then sorts these lists and calculates the difference
 * between occurances of the k-mers in A and B.
 */
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

