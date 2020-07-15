#include <random>

#include <benchmark/benchmark.h>

#include <seqan/sequence.h>
#include <seqan/index.h>

template <typename alphabet_t>
seqan::String<alphabet_t> generate_sequence_seqan2(size_t const len = 500,
                                                   size_t const variance = 0,
                                                   size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, seqan::ValueSize<alphabet_t>::VALUE - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    seqan::String<alphabet_t> sequence;
    size_t length = dis_length(gen);
    seqan::resize(sequence, length);

    for (size_t l = 0; l < length; ++l)
        sequence[l] = alphabet_t{dis_alpha(gen)};

    return sequence;
}

template <typename TAlphabet>
void rank(benchmark::State & state)
{
    size_t const sigma = seqan::ValueSize<TAlphabet>::VALUE;
    size_t const text_size{static_cast<size_t>(state.range(0))};

    if (state.range(0) <= static_cast<decltype(state.range(0))>(0))
        throw std::invalid_argument{"Text needs to be at least 1 character long."};

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<> sigma_dist(0, sigma - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<size_t> position_dist(0, text_size - 2);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    seqan::String<TAlphabet> const text = generate_sequence_seqan2<TAlphabet>(text_size);

    using epr_config = seqan::LevelsPrefixRDConfig<size_t, seqan::Alloc<>, 2, 1>;
    using level_config = seqan::Levels<void, epr_config>;
    seqan::RankDictionary<TAlphabet, level_config> epr(text);

    for (auto _ : state)
    {
        size_t const pos{position_gen()};
        TAlphabet const val{sigma_gen()};

        benchmark::DoNotOptimize(getRank(epr, pos, val));
    }
}


BENCHMARK_TEMPLATE(rank, seqan::Dna)->RangeMultiplier(100)->Range(100, 100'000'000);
BENCHMARK_TEMPLATE(rank, seqan::AminoAcid)->RangeMultiplier(100)->Range(100, 100'000'000);

BENCHMARK_MAIN();


