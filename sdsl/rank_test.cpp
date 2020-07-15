#include <gtest/gtest.h>
#include <sdsl/rank_support_int.hpp>

#include <random>

TEST(rank, crash)
{
    constexpr uint8_t sigma{4u};

    std::mt19937_64 engine(12345);
    std::uniform_int_distribution<> dist(0, sigma - 1);
    auto gen = [&dist, &engine]() { return dist(engine); };

    sdsl::int_vector<> text(10'000, 0, 2);
    std::generate(text.begin(), text.end(), gen);

    std::vector<size_t> cnt_prefix_rank(sigma, 0);
    sdsl::rank_support_int_v<sigma> rank_vector(&text);

    rank_vector.rank(10'000, 3);
}

TEST(rank, rankier)
{
    constexpr uint8_t sigma{32u};

    std::mt19937_64 engine(12345);
    std::uniform_int_distribution<> dist(0, sigma - 1);
    auto gen = [&dist, &engine]() { return dist(engine); };

    sdsl::int_vector<> text(500, 0, 5);
    std::generate(text.begin(), text.end(), gen);

    std::vector<size_t> cnt_prefix_rank(sigma, 0);
    sdsl::rank_support_int_v<sigma> rank_vector(&text);

    for (size_t i = 0; i <= text.size(); ++i)
    {
        for (size_t c = 0; c < sigma; ++c)
        {
            if (i > 0 && text[i - 1] <= c)
                ++cnt_prefix_rank[c];

            if (c > 0)
                ASSERT_EQ(cnt_prefix_rank[c] - cnt_prefix_rank[c - 1], rank_vector.rank(i, c)) << "i=" << i << " c=" << c << "/" << sigma - 1;
            else
                ASSERT_EQ(cnt_prefix_rank[c], rank_vector.rank(i, c)) << "i=" << i << " c=" << c << "/" << sigma - 1;
        }
    }
}
