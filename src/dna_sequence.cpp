#include "dna_sequence.hpp"

namespace DNA {
	Sequence::Sequence(std::string const& dna_data)
		: _size{dna_data.size()}
		, _entries(new uint8_t[(dna_data.size()+3)>>2]())
	{
		for (size_t i {}; i < dna_data.size(); ++i)
			_entries[i >> 2] |= to_value(dna_data[i]) << (2 * (i & 3));
	}

	Sequence::~Sequence() { delete[] _entries; }

	std::ostream& operator<<(std::ostream& os, Sequence const& sequence) 
	{
		for (size_t i {}; i < sequence.size(); ++i)
			os << Sequence::to_letter(sequence[i]);

		return os;
	}
}
