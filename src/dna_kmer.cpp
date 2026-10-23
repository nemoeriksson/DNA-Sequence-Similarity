#include "dna_kmer.hpp"
#include "dna_sequence.hpp"

namespace DNA
{
KmerArray::KmerArray(Sequence const& seq_a,
                     Sequence const& seq_b,
                     const size_t k)
    : _differences(new difference_type[1ULL << (2 * k)]())
{
    uint64_t kmer_a{};
    uint64_t kmer_b{};

    for (size_t i{}; i < k - 1; ++i)
    {
        kmer_a = (kmer_a << 2) | seq_a[i];
        kmer_b = (kmer_b << 2) | seq_b[i];
    }

    for (size_t i{ k - 1 }; i < seq_a.size(); ++i)
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
}
