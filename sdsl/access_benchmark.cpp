#include <benchmark/benchmark.h>
#include <sdsl/wt_epr.hpp>

#include <algorithm>
#include <cmath>
#include <random>

template <uint8_t sigma>
void access_block(benchmark::State & state)
{
    if (state.range(0) <= static_cast<decltype(state.range(0))>(0))
        throw std::invalid_argument{"Text needs to be at least 1 character long."};

    uint8_t const log_sigma{static_cast<uint8_t>(std::clamp(std::ceil(std::log2(sigma)), 1.0, 64.0))};
    size_t const text_size{static_cast<size_t>(state.range(0))};

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<> sigma_dist(0, sigma - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    sdsl::epr::int_vector<> text(text_size, 0, log_sigma);
    std::generate(text.begin(), text.end(), sigma_gen);
    sdsl::epr::rank_support_int_v<sigma> rank_support(&text);

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<> position_dist(0, rank_support.m_block.size() - 1);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(rank_support.m_block[position_gen()]);
    }
}

template <uint8_t sigma>
void access_superblock(benchmark::State & state)
{
    if (state.range(0) <= static_cast<decltype(state.range(0))>(0))
        throw std::invalid_argument{"Text needs to be at least 1 character long."};

    uint8_t const log_sigma{static_cast<uint8_t>(std::clamp(std::ceil(std::log2(sigma)), 1.0, 64.0))};
    size_t const text_size{static_cast<size_t>(state.range(0))};

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<> sigma_dist(0, sigma - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    sdsl::epr::int_vector<> text(text_size, 0, log_sigma);
    std::generate(text.begin(), text.end(), sigma_gen);
    sdsl::epr::rank_support_int_v<sigma> rank_support(&text);

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<> position_dist(0, rank_support.m_superblock.size() - 1);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(rank_support.m_superblock[position_gen()]);
    }
}

// BENCHMARK_TEMPLATE(access_block, 4)->RangeMultiplier(100)->Range(100, 1'000'000'000);
// BENCHMARK_TEMPLATE(access_block, 27)->RangeMultiplier(100)->Range(100, 1'000'000'000);
BENCHMARK_TEMPLATE(access_superblock, 4)->RangeMultiplier(100)->Range(100, 1'000'000'000);
BENCHMARK_TEMPLATE(access_superblock, 27)->RangeMultiplier(100)->Range(100, 1'000'000'000);
// BENCHMARK_TEMPLATE(rank, 8)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 16)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 32)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 64)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 128)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 255)->RangeMultiplier(100)->Range(100, 1'000'000);

BENCHMARK_MAIN();
