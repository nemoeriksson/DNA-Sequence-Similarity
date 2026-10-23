#include "dna_similarity.hpp"
#include "dna_kmer.hpp"
#include "dna_sequence.hpp"

/* Private declarations */

/*
 * Uses an array with one entry per possible k-mer (length = 4^k)
 * where the n:th entry is equal to abs(C_A(x) - C_B(x)) where
 * C_A(x) The count of occurrences of the k-mer x inside DNA sequence A
 * C_B(x) The count of occurrences of the k-mer x inside DNA sequence B
 *
 * Recommended for lower values k since requierd memory quickly increases
 */
static double calculate_using_kmer_array(DNA::Sequence const& seq_a,
                                         DNA::Sequence const& seq_b,
                                         const size_t k);

/*
 * Adds all k-mers into a sorted list with its occurrences in separate lists.
 */
template<typename T>
static double calculate_using_kmer_sorted(DNA::Sequence const& seq_a,
                                          DNA::Sequence const& seq_b,
                                          const size_t k);

template<typename T>
static size_t get_count(DNA::KmerSorted<T> const& kmer_sorted,
                        T target,
                        size_t index,
                        T (DNA::KmerSorted<T>::*at)(size_t) const);

/* Public definitions */

namespace DNA
{
double get_similarity(Sequence const& seq_a,
                      Sequence const& seq_b,
                      const size_t k)
{
    if (k <= 8)
        return calculate_using_kmer_array(seq_a, seq_b, k);

    else if (k <= 16)
        return calculate_using_kmer_sorted<uint32_t>(seq_a, seq_b, k);

    return calculate_using_kmer_sorted<uint64_t>(seq_a, seq_b, k);
}
}

/* Private definitions */

// TODO: Add vectorization where applicable
static double calculate_using_kmer_array(DNA::Sequence const& seq_a,
                                         DNA::Sequence const& seq_b,
                                         const size_t k)
{
    const DNA::KmerArray kmer_arr{ seq_a, seq_b, k };
    const size_t max_kmer{ 1ULL << (2 * k) };
    const size_t mid = seq_a.size() + seq_b.size() - 2 * k + 2;

    uint64_t diff{ 0 };
    for (size_t i{}; i < max_kmer; ++i)
        diff += std::abs(kmer_arr[i]);

    return (double)(mid - diff) / (double)(mid + diff);
}

template<typename T>
static double calculate_using_kmer_sorted(DNA::Sequence const& seq_a,
                                          DNA::Sequence const& seq_b,
                                          const size_t k)
{
    const DNA::KmerSorted<T> kmer_sorted{ seq_a, seq_b, k };
    const size_t mid = seq_a.size() + seq_b.size() - 2 * k + 2;

    size_t diff{ 0 };
    size_t ai{ 0 };
    size_t bi{ 0 };

    while (ai < kmer_sorted.size() && bi < kmer_sorted.size())
    {
        T kmer_a{ kmer_sorted.at_a(ai) };
        T kmer_b{ kmer_sorted.at_b(bi) };

        if (kmer_a < kmer_b)
        {
            size_t count =
                get_count(kmer_sorted, kmer_a, ai, &DNA::KmerSorted<T>::at_a);
            diff += count;
            ai += count;
        }
        else if (kmer_a > kmer_b)
        {
            size_t count =
                get_count(kmer_sorted, kmer_b, bi, &DNA::KmerSorted<T>::at_b);
            diff += count;
            bi += count;
        }
        else
        {
            size_t count_a =
                get_count(kmer_sorted, kmer_a, ai, &DNA::KmerSorted<T>::at_a);
            size_t count_b =
                get_count(kmer_sorted, kmer_b, bi, &DNA::KmerSorted<T>::at_b);
            diff += std::abs((int)count_a - (int)count_b);
            ai += count_a;
            bi += count_b;
        }
    }

    // Add those left over
    while (ai < kmer_sorted.size())
    {
        T kmer = kmer_sorted.at_a(ai);
        size_t count =
            get_count(kmer_sorted, kmer, ai, &DNA::KmerSorted<T>::at_a);
        diff += count;
        ai += count;
    }
    while (bi < kmer_sorted.size())
    {
        T kmer = kmer_sorted.at_b(bi);
        size_t count =
            get_count(kmer_sorted, kmer, bi, &DNA::KmerSorted<T>::at_b);
        diff += count;
        bi += count;
    }

    return (double)(mid - diff) / (double)(mid + diff);
}
template<typename T>
static size_t get_count(DNA::KmerSorted<T> const& kmer_sorted,
                        T target,
                        size_t index,
                        T (DNA::KmerSorted<T>::*at)(size_t) const)
{
    size_t count{};
    while (index + count < kmer_sorted.size() &&
           (kmer_sorted.*at)(index + count) == target)
    {
        ++count;
    }

    return count;
}
