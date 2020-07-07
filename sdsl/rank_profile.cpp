#include <sdsl/wt_epr.hpp>

#include <algorithm>
#include <cmath>
#include <random>

template <uint8_t sigma>
size_t rank(size_t const text_size, size_t const num_executions)
{
    uint8_t const log_sigma{static_cast<uint8_t>(std::clamp(std::ceil(std::log2(sigma)), 1.0, 64.0))};
    // size_t const text_size{static_cast<size_t>(state.range(0))};

    using size_type = typename sdsl::epr::rank_support_int_v<sigma>::size_type;
    using value_type = typename sdsl::epr::rank_support_int_v<sigma>::value_type;

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<value_type> sigma_dist(0, sigma - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    sdsl::epr::int_vector<log_sigma> text(text_size, 0);
    std::generate(text.begin(), text.end(), sigma_gen);
    sdsl::epr::rank_support_int_v<sigma,log_sigma> rank_support(&text);

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<size_type> position_dist(0, text_size - 1);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    size_t rank = 0;
    for (size_t i = 0; i < num_executions; ++i)
    {
        size_type const position{position_gen()};
        value_type const value{sigma_gen()};
        rank += rank_support.rank(position, value);
    }

    return rank;
}

int main()
{
    size_t r = rank<4>(1'000'000, 1'000'000'000);
    std::cout << "rank: " << r << "\n";
}
