#include <benchmark/benchmark.h>
#include <sdsl/wt_epr.hpp>

#include <algorithm>
#include <cmath>
#include <random>

template <uint8_t sigma>
void rank(benchmark::State & state)
{
    if (state.range(0) <= static_cast<decltype(state.range(0))>(0))
        throw std::invalid_argument{"Text needs to be at least 1 character long."};

    uint8_t const log_sigma{static_cast<uint8_t>(std::clamp(std::ceil(std::log2(sigma)), 1.0, 64.0))};
    size_t const text_size{static_cast<size_t>(state.range(0))};

    using size_type = typename sdsl::rank_support_int_v<sigma>::size_type;
    using value_type = typename sdsl::rank_support_int_v<sigma>::value_type;

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<value_type> sigma_dist(0, sigma - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    sdsl::int_vector<> text(text_size, 0);
    std::generate(text.begin(), text.end(), sigma_gen);
    sdsl::rank_support_int_v<sigma> rank_support(&text);

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<size_type> position_dist(0, text_size - 1);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    for (auto _ : state)
    {
        size_type const position{position_gen()};
        value_type const value{sigma_gen()};
        benchmark::DoNotOptimize(rank_support.rank(position, value));
    }
}

BENCHMARK_TEMPLATE(rank, 4)->RangeMultiplier(100)->Range(100, 100'000'000);
// BENCHMARK_TEMPLATE(rank, 4)->Arg(100000000);
BENCHMARK_TEMPLATE(rank, 27)->RangeMultiplier(100)->Range(100, 100'000'000);
// BENCHMARK_TEMPLATE(rank, 8)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 16)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 32)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 64)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 128)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, 255)->RangeMultiplier(100)->Range(100, 1'000'000);

BENCHMARK_MAIN();
