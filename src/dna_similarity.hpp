#ifndef DNA_SIM_HPP
#define DNA_SIM_HPP

#include "dna_sequence.hpp"
#include "iostream"

namespace DNA
{
/*
 * Calculates the similarity between two DNA sequences sing k-map overlap
 * and the Jaccard index.
 */
double get_similarity(Sequence const& seq_a,
                      Sequence const& seq_b,
                      const size_t k);
}

#endif
