#include "dna_similarity.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

static constexpr std::string generate_dna_sequence(size_t length);
static inline void mutate_dna_sequence(std::string& dna_str,
                                       const size_t mutations);

int main()
{
    constexpr size_t min_dna_length_exp{ 12 };
    constexpr size_t max_dna_length_exp{ 22 };
    constexpr size_t max_k{ 32 };
    constexpr double zero{ 1e-8 };

    // Output CSV header
    std::cout << "Sequence Length;"
              << "Mutations;"
              << "k;"
              << "Similarity Score;"
              << "Avg Time;" << std::endl;
    // Warmup
    double a{ 1 };
    for (size_t i{}; i < 100000000; ++i)
        a += rand();
    std::cerr << "Warmup complete " << a << std::endl;

    for (size_t exp{ min_dna_length_exp }; exp <= max_dna_length_exp; ++exp)
    {
        const size_t sequence_length = 1 << exp;
        const size_t mutations_per_iteration{ sequence_length / 100 };
        const std::string dna_sequence_str{ generate_dna_sequence(
            sequence_length) };
        std::string mutated_dna_string{ dna_sequence_str };

        double similarity{};
        size_t mutations{};
        do
        {
            for (size_t k{ 1 }; k <= max_k; ++k)
            {
                auto start_time = std::chrono::system_clock::now();

                similarity = DNA::get_similarity(
                    dna_sequence_str, mutated_dna_string, k);

                auto end_time = std::chrono::system_clock::now();
                auto elapsed = end_time - start_time;
                double avg_ms =
                    std::chrono::duration<double, std::milli>(elapsed).count();

                // Output CSV data
                std::cout << std::format("{:9};{:9};{:9};{:.18f};{:.8f}\n",
                                         sequence_length,
                                         mutations,
                                         k,
                                         similarity,
                                         avg_ms);
            }
            mutate_dna_sequence(mutated_dna_string, mutations_per_iteration);
            mutations += mutations_per_iteration;
        } while (similarity > zero);
    }

    return 0;
}

static constexpr std::string generate_dna_sequence(const size_t length)
{
    std::string dna_str;

    for (size_t i{}; i < length; ++i)
        dna_str += DNA::Sequence::to_letter(std::rand() % 4);

    return dna_str;
}

static inline void mutate_dna_sequence(std::string& dna_str,
                                       const size_t mutations)
{
    size_t size{ dna_str.size() };
    for (size_t i{}; i < mutations; ++i)
    {
        size_t index = rand() % size;
        dna_str[index] += 1 << 1;
    }
}
