#include "dna_kmer.hpp"
#include "dna_sequence.hpp"
#include <algorithm>

static void add_to_bucket(
    std::vector<DNA::KmerBuckets::BucketItem>& bucket,
    DNA::KmerBuckets::BucketItem const&& item);

namespace DNA
{
KmerArray::KmerArray(Sequence const& seq_a,
                     Sequence const& seq_b,
                     const size_t k)
    : _differences(new difference_type[1ULL << (2 * k)]())
{
    uint64_t kmer_a{};
    uint64_t kmer_b{};

    for (size_t i{}; i < k-1; ++i)
    {
        kmer_a = (kmer_a << 2) | seq_a[i];
        kmer_b = (kmer_b << 2) | seq_b[i];
    }

    for (size_t i{ k-1 }; i < seq_a.size(); ++i)
    {
        kmer_a = ((kmer_a << 2) | seq_a[i]) & ((1 << (2 * k)) - 1);
        kmer_b = ((kmer_b << 2) | seq_b[i]) & ((1 << (2 * k)) - 1);

        ++_differences[kmer_a];
        --_differences[kmer_b];
    }
}

KmerArray::~KmerArray()
{
    delete[] _differences;
}

KmerBuckets::KmerBuckets(Sequence const& seq_a,
                         Sequence const& seq_b,
                         const size_t k)

{
    for (size_t i{}; i < 256; ++i)
        _buckets[i].reserve(4096);

    uint64_t kmer_a{};
    uint64_t kmer_b{};

    for (size_t i{}; i < k-1; ++i)
    {
        kmer_a = (kmer_a << 2) | seq_a[i];
        kmer_b = (kmer_b << 2) | seq_b[i];
    }

    for (size_t i{ k-1 }; i < seq_a.size(); ++i)
    {
        kmer_a = ((kmer_a << 2) | seq_a[i]);
        kmer_b = ((kmer_b << 2) | seq_b[i]);

        if (k != 32) // If k=32 the entire kmer fills the 64
                     // bits so no mask is needed to extract it
        {
            kmer_a &= ((1ULL << (2 * k)) - 1);
            kmer_b &= ((1ULL << (2 * k)) - 1);
        }

        size_t bucket_index_a = kmer_a >> (2 * k - 8);
        size_t bucket_index_b = kmer_b >> (2 * k - 8);
        uint64_t kmer_a_seg = kmer_a & ((1ULL << (2 * k - 8)) - 1);
        uint64_t kmer_b_seg = kmer_b & ((1ULL << (2 * k - 8)) - 1);

        add_to_bucket(_buckets[bucket_index_a], BucketItem{0, 1, kmer_b_seg});
        add_to_bucket(_buckets[bucket_index_b], BucketItem{1, 0, kmer_a_seg});

        //_buckets[bucket_index_a].emplace_back(1, 0, kmer_a_seg);
        //_buckets[bucket_index_b].emplace_back(0, 1, kmer_b_seg);
    }

	/*
    for (size_t i{}; i < 256; ++i)
        std::sort(_buckets[i].begin(),
                  _buckets[i].end(),
                  [](BucketItem const& a, BucketItem const& b)
                  { return a.kmer_segment < b.kmer_segment; });
				  */
}
}

static void add_to_bucket(
    std::vector<DNA::KmerBuckets::BucketItem>& bucket,
    DNA::KmerBuckets::BucketItem const&& item)
{
	auto insert_index = std::lower_bound(
			bucket.begin(),
			bucket.end(),
			item,
			[](auto const& a, auto const& b) {
				return a.kmer_segment < b.kmer_segment;});

	if (insert_index != bucket.begin())
	{
		auto& previous {*(insert_index-1)};
		
		if (previous.kmer_segment == item.kmer_segment)
		{
			if (previous.count_a + item.count_a <= 16
					&& previous.count_b + item.count_b <= 16)
			{
				previous.count_a += item.count_a;
				previous.count_b += item.count_b;
			}
			return;
		}
	}

	bucket.insert(insert_index, item);
}
