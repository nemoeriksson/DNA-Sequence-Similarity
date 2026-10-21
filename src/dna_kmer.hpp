#ifndef DNA_KMER_HPP
#define DNA_KMER_HPP

#include "dna_sequence.hpp"
#include <cstring>

namespace DNA
{
using difference_type = int16_t;

class KmerArray
{
public:
    KmerArray(Sequence const& seq_a, Sequence const& seq_b, const size_t k);

    ~KmerArray();

    inline difference_type operator[](const size_t i) const
    {
        return _differences[i];
    }

private:
    difference_type* const _differences;
};

template<typename T> class KmerSorted
{
public:
    KmerSorted(Sequence const& seq_a, Sequence const& seq_b, const size_t k)
			: _size(seq_a.size() - k + 1)
			, _kmers_a(new T[_size]())
			, _kmers_b(new T[_size]())
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
            kmer_a = ((kmer_a << 2) | seq_a[i]);
            kmer_b = ((kmer_b << 2) | seq_b[i]);

			if (k != 32) // If k=32 the entire kmer fills the 64
                     // bits so no mask is needed to extract it
			{
				kmer_a &= ((1ULL << (2 * k)) - 1);
				kmer_b &= ((1ULL << (2 * k)) - 1);
			}

			_kmers_a[i-k+1] = kmer_a;
			_kmers_b[i-k+1] = kmer_b;
        }

		std::sort(_kmers_a, _kmers_a+_size);
		std::sort(_kmers_b, _kmers_b+_size);
    }

    ~KmerSorted()
    {
        delete[] _kmers_a;
        delete[] _kmers_b;
    }

    inline size_t size() const { return _size; }

    inline T at_a(const size_t i) const { return _kmers_a[i]; }
    inline T at_b(const size_t i) const { return _kmers_b[i]; }

private:
	size_t _size;
	T* _kmers_a;
	T* _kmers_b;
};
}

#endif
