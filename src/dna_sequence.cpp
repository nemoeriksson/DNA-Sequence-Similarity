#include "dna_sequence.hpp"

namespace DNA {
	Sequence::~Sequence()
	{
		delete[] _entries;
	}

	std::ostream& operator<<(std::ostream& os, Sequence const& sequence) 
	{
		for (size_t i {}; i < sequence.size(); ++i)
		{
			os << Sequence::to_letter(sequence[i]);
		}

		return os;
	}
}
