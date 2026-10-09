#ifndef DNA_SIM_HPP
#define DNA_SIM_HPP

#include "iostream"
#include "dna_sequence.hpp"

namespace DNA
{
	double get_similarity(Sequence const& seq_a, 
			Sequence const &seq_b,
			const size_t k);
}

#endif
