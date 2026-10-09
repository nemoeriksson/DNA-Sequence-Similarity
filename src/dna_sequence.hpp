#ifndef DNA_SEQ
#define DNA_SEQ

#include <iostream>
#include <cstdint>
#include <ostream>
#include <string>

namespace DNA
{
	class Sequence
	{
		public:
			constexpr Sequence(std::string const& dna_data)
				: _size{dna_data.size()}
				, _entries(new uint8_t[(dna_data.size()+3)>>2]())
			{
				for (size_t i {}; i < dna_data.size(); ++i)
					_entries[i >> 2] |= to_value(dna_data[i]) << (2 * (i & 3));
			}

			~Sequence();

			inline constexpr size_t size() const { return _size; } 
		
			inline constexpr uint8_t operator[](const size_t i) const 
			{
				return (_entries[i >> 2] >> (2 * (i & 3))) & 0b11;
			}

			inline constexpr size_t get_kmer_count(const size_t k) const
			{
				return size()-k+1;
			}

			static inline constexpr char to_letter(const uint8_t entry)
			{
				static constexpr char map[4] {'A', 'C', 'T', 'G'};
				return map[entry&0b11];
			}

			static inline constexpr uint8_t to_value(const char letter)
			{
				return (letter>>1)&0b11;
			}
		private:
			const size_t _size {};
			uint8_t * const _entries;

			// Size: 16B
			// Alignment: 8B
	};

	std::ostream& operator<<(std::ostream& os, Sequence const& sequence);
}

#endif
