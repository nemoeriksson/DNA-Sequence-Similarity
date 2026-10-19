#ifndef DNA_SEQ_HPP
#define DNA_SEQ_HPP

#include <iostream>

namespace DNA
{
/*
 * Stores a sequence of DNA bases (A, C, T, or G) bit-packed to take
 * 2 bits each.
 *
 * ===== Stats & Storage Usage (*_entries) =====
 * Size = ceil(n / 256) bytes
 *
 * Fills one cache line (64B) at size 256.
 *
 * My device:
 * Fills 100% of L1 cache at n =   131072 (32 KiB)
 *		 100%    L2			 n =  4194304 ( 1 MiB)
 *		 100%	 L3			 n = 67108864 (16 MiB)
 *
 * University device (salt):
 * Fills 100% of L1 cache at n =   131072 (32 KiB)
 *		 100%    L2			 n =  1048576 (256 KiB)
 *		 100%	 L3			 n = 67108864 ( 16 MiB)
 *
 */
class Sequence
{
public:
    Sequence(std::string const& dna_data);

    ~Sequence();

    inline size_t size() const
    {
        return _size;
    }

    inline uint8_t operator[](const size_t i) const
    {
        return (_entries[i >> 2] >> (2 * (i & 3))) & 0b11;
    }

    inline size_t get_kmer_count(const size_t k) const
    {
        return size() - k + 1;
    }

    inline uint8_t* get_entries() const
    {
        return _entries;
    }

    static inline constexpr char to_letter(const uint8_t entry)
    {
        static constexpr char map[4]{ 'A', 'C', 'T', 'G' };
        return map[entry & 0b11];
    }

    static inline constexpr uint8_t to_value(const char letter)
    {
        return (letter >> 1) & 0b11;
    }

private:
    const size_t _size{};
    uint8_t* const _entries;
};

std::ostream& operator<<(std::ostream& os, Sequence const& sequence);
}

#endif
