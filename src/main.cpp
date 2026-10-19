#include "dna_sequence.hpp"
#include "dna_similarity.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>

static constexpr std::string generate_dna_sequence(size_t length);

int main()
{
    constexpr size_t dna_sequence_length{ 1 << 20 };
    constexpr size_t iterations{ 1 };
    constexpr size_t max_k{ 32 };

    DNA::Sequence seq_a{ generate_dna_sequence(dna_sequence_length) };
    DNA::Sequence seq_b{ generate_dna_sequence(dna_sequence_length) };

    // Output CSV header
    std::cout << "Sequence Length;"
              << "k;"
              << "Similarity Score;"
              << "Average Calculation Time;" << std::endl;

    double similarity;
    for (size_t k{ 1 }; k <= max_k; ++k)
    {
        auto start_time = std::chrono::system_clock::now();
        for (size_t i{}; i < iterations; ++i)
            similarity = DNA::get_similarity(seq_a, seq_b, k);
        auto end_time = std::chrono::system_clock::now();

        auto elapsed = end_time - start_time;
        double avg_ms =
            std::chrono::duration<double, std::milli>(elapsed / iterations)
                .count();

        // CSV data
        std::cout << dna_sequence_length << ';';
        std::cout << k << ';';
        std::cout << similarity << ';';
        std::cout << std::format("{:.2f}ms", avg_ms) << std::endl;
    }

    return 0;
}

static constexpr std::string generate_dna_sequence(size_t length)
{
    std::string dna_str;

    for (size_t i{}; i < length; ++i)
        dna_str += DNA::Sequence::to_letter(std::rand() % 4);

    return dna_str;
}
