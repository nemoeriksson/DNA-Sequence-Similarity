#ifndef DNA_SIM_HPP
#define DNA_SIM_HPP

#include "dna_sequence.hpp"
#include "iostream"

namespace DNA
{
/*
 * Calculates the similarity between two DNA sequences using k-map overlap
 * and the Jaccard index.
 */
double get_similarity(std::string const& dna_seq_str_a,
                      std::string const& dna_seq_str_b,
                      const size_t k);
}

#endif
