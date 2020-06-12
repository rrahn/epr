// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int_v.hpp
    \brief rank_support_int_v.hpp contains rank_support_int_v.
	\author Christopher Pockrandt
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT_V
#define INCLUDED_SDSL_RANK_SUPPORT_INT_V

#include "rank_support_int.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl {
namespace epr {

//! A rank structure proposed by Christopher Pockrandt
/*!
 * This data structure is similar to rank data structures on bit vectors.
 * It supports constant time rank and prefix_rank queries on int vectors.
 *
 * \tparam alphabet_size         Size of the alphabet represented in the int_vector, i.e., largest value + 1.
 * \tparam words_per_block       Words per block (equivalent to the number of popcount operations in the worst-case per rank query).
 * \tparam blocks_per_superblock Blocks per superblock.
 *
 * \par Reference
 *    Christopher Pockrandt:
 *    EPR-Dictionaries: A practical and fast data structure for constant time searches in unidirectional and bidirectional FM-indices.
 *    WEA 2008: 154-168
 *
 * @ingroup rank_support_group
 */
template <uint8_t alphabet_size, uint8_t vector_width=0, uint8_t words_per_block = 1, uint8_t blocks_per_superblock = 4>
class rank_support_int_v : public rank_support_int<alphabet_size,vector_width> {
public:
	// typedef int_vector<vector_width> int_vector_type;
	typedef typename rank_support_int<alphabet_size,vector_width>::size_type size_type;
	typedef typename rank_support_int<alphabet_size,vector_width>::value_type value_type;

public:
	int_vector<8> m_block;
	int_vector<32> m_superblock;

	static constexpr uint64_t values_per_word{64ULL / rank_support_int<alphabet_size,vector_width>::sigma_bits};
	static constexpr uint32_t values_per_block{words_per_block * values_per_word};

public:
	explicit rank_support_int_v(const int_vector<vector_width>* v = nullptr) : rank_support_int<alphabet_size,vector_width>(v)
	{
	    static_assert(blocks_per_superblock > 1, "There must be at least two blocks per superblock!");
		constexpr uint8_t max_letter{static_cast<uint8_t>(this->sigma) - 1};

		if (v == nullptr)
		{
			return;
		}
		else if (v->empty())
		{
			m_block.resize(max_letter, 0);
			m_superblock.resize(max_letter, 0);
			return;
		}

		constexpr uint64_t words_per_superblock{words_per_block * blocks_per_superblock};
		constexpr uint64_t values_per_superblock{blocks_per_superblock * values_per_block};
		constexpr uint64_t new_width{ceil_log2(values_per_superblock)};
		// std::cout << new_width << '\n';
		m_block.width(new_width);
		// m_superblock.width(64);

		// NOTE: number of elements is artificially increased by one because rank can be called on m_v[size()]
		uint64_t const word_count = ((this->m_v->size() - 1 + 1) / values_per_word) + 1; // equivalent to ceil(m_v->size() / values_per_word)
		uint64_t const block_count = ((word_count - 1) / words_per_block) + 1; // equivalent to ceil(word_count / words_per_block)

		// for each superblock we only need `blocks_per_superblock-1` instead of `blocks_per_superblock` blocks.
		// for the last superblock we can subtract the last unused blocks.
        size_type const blocks_needed = (((block_count - 1) / blocks_per_superblock) + 1) * (blocks_per_superblock - 1)
									  - ((blocks_per_superblock - (block_count % blocks_per_superblock)) % blocks_per_superblock);
		size_type const block_size = (blocks_needed + 1) * max_letter;
		size_type const superblock_size = (((word_count - 1) / words_per_superblock) + 1) * max_letter; // equivalent to ceil(word_count / words_per_superblock) * max_letter
		m_block.resize(block_size);
		m_superblock.resize(superblock_size);

		uint64_t const * data = this->m_v->data();
		std::vector<uint64_t> buf_blocks(max_letter, 0);
		std::vector<uint64_t> buf_superblocks(max_letter, 0);

		for (uint64_t v = 0; v < max_letter; ++v)
			m_superblock[v] = 0;

		// Precompute blocks and superblocks
		// NOTE: divisors in modulo operations are constexpr and hence are expected to be cheap
		for (uint64_t word_id = 0, block_id = 0, superblock_id = max_letter; word_id < word_count; ++word_id)
		{
			for (uint64_t v = 0; v < max_letter; ++v)
				buf_blocks[v] += this->full_word_prefix_rank(data, word_id, v);

			// counted the values in the last word of the current block
			if (word_id % words_per_block == (words_per_block - 1))
			{
                if (word_id % words_per_superblock != (words_per_superblock - 1))
                {
                    if (block_id < m_block.size()) // TODO: rewrite for loop to eliminate need for if clause here
                    {
                        for (uint64_t v = 0; v < max_letter; ++v)
                            m_block[block_id + v] = buf_blocks[v];
                        block_id += max_letter;
                    }
                }
                else
                { // don't store block information for the last block in the superblock!
                    if (superblock_id < m_superblock.size()) // TODO: rewrite for loop to eliminate need for if clause here
                    {
                        for (uint64_t v = 0; v < max_letter; ++v)
                        {
                            buf_superblocks[v] += buf_blocks[v];
                            m_superblock[superblock_id + v] = buf_superblocks[v];
                            buf_blocks[v] = 0; // reset blocks
                        }
                    }
                    superblock_id += max_letter;
                }
			}
		}

        // std::cout << "\nBlocks:\n";
        // for (uint64_t i = 0; i < m_block.size(); i += max_letter)
		// {
		// 	for (uint64_t v = 0; v < max_letter; ++v)
	    //         std::cout << (unsigned)m_block[i + v] << ' ';
        //     std::cout << "| ";
		// }
        // std::cout << "\nSuperBlocks:\n";
        // for (uint64_t i = 0; i < m_superblock.size(); i += max_letter)
		// {
		// 	for (uint64_t v = 0; v < max_letter; ++v)
        //     	std::cout << (unsigned)m_superblock[i + v] << ' ';
        //     std::cout << "| ";
		// }
        // std::cout << "\n\n";
	}

	rank_support_int_v(const rank_support_int_v&) = default;
	rank_support_int_v(rank_support_int_v&&)	  = default;
	rank_support_int_v& operator=(const rank_support_int_v&) = default;
	rank_support_int_v& operator=(rank_support_int_v&&) = default;

	//! Counts the occurrences of v in the prefix [0..idx-1]
	/*! \param idx Argument for the length of the prefix v[0..idx-1].
		 *  \param v Argument which value to count.
	     *  \sa prefix_rank
	     */
	size_type rank(const size_type idx, const value_type v) const
	{
		switch (v)
		{
			case 0 : return prefix_rank(idx, v);
			case this->sigma-1 : return idx - prefix_rank(idx, v - 1);
			// case this->sigma-1 : return prefix_rank(idx, v - 1);
			default: return prefix_rank2(idx, v);
			// default: return prefix_rank(idx, v);
		}
	}

	//! Alias for rank(idx, v)
	inline size_type operator()(const size_type idx, const value_type v) const { return rank(idx, v); }

	//! Counts the occurrences of elements smaller or equal to v in the prefix [0..idx-1]
	/*! \param idx Argument for the length of the prefix v[0..idx-1].
		 *  \param v Argument which value (including smaller values) to count.
	     *  \sa rank
	     */
	size_type prefix_rank(const size_type idx, const value_type v) const
	{
		assert(this->m_v != nullptr);
		assert(idx <= this->m_v->size());
        assert(v <= this->sigma);

		if (unlikely(v == this->sigma - 1)) // TODO actually
		{
			return idx;
		}

		constexpr uint8_t max_letter{this->sigma - 1};

		size_type const block_id{idx / values_per_block};
		size_type const superblock_id{block_id / blocks_per_superblock};
        size_type const block_id_in_superblock{block_id - blocks_per_superblock * superblock_id};
      //size_type const block_id_in_superblock{block_id % blocks_per_superblock};

		// retrieve superblock value
        size_type res = m_superblock[max_letter * superblock_id + v];
		// std::cout << "res1=" << res << '\n';

		// retrieve block value
		size_type const cache = block_id_in_superblock > 0;
		size_type const pppp = max_letter * (superblock_id * (blocks_per_superblock - 1) + (block_id_in_superblock + !cache - 1)) + v;
		res += cache * m_block[pppp];

        // if (block_id_in_superblock > 0)
            // res += m_block[max_letter * superblock_id * (blocks_per_superblock - 1) + (block_id_in_superblock - 1)  * max_letter + v];

		// std::cout << "res2=" << res << '\n';
		// compute in-block queries for all words before the in-block queries
		// this only applies when multiple words are in one block
        if constexpr (words_per_block > 1)
        {
            size_type const word_id{idx / values_per_word};
            uint64_t w{word_id - (word_id % words_per_block)};
            while (w < word_id)
            {
                res += this->full_word_prefix_rank(this->m_v->data(), w, v);
                ++w;
            }
			// std::cout << "res3=" << res << '\n';
        }

		// compute in-block query
		// std::cout << "values_per_block=" << values_per_block << '\n';
		// if (idx % values_per_block != 0)
		// {
		size_type cache2 = idx % values_per_block != 0;
		// size_type const cache2 = 1;
		res += cache2 * this->word_prefix_rank(this->m_v->data(), idx, v);
			// std::cout << "idx=" << idx << '\n';
		// }
		// std::cout << "res4=" << res << '\n';

		return res;
	}

	size_type prefix_rank2(const size_type idx, const value_type v) const
	{
		assert(this->m_v != nullptr);
		assert(idx <= this->m_v->size());
        assert(v <= this->sigma);

		constexpr uint8_t const max_letter{this->sigma - 1};

		size_type const block_id{idx / values_per_block};
		size_type const superblock_id{block_id / blocks_per_superblock};
		size_type const block_id_in_superblock{block_id - blocks_per_superblock * superblock_id};

		size_type const access1 = max_letter * superblock_id + v;
		auto it1 = m_superblock.begin() + (access1 - 1);
		size_type res_lower = *it1;
		++it1;
		size_type res_upper = *it1;
        // size_type res_upper = m_superblock[access1];
        // size_type res_lower = m_superblock[access1 - 1];

		size_type const cache = block_id_in_superblock > 0;
		size_type const access2 = max_letter * (superblock_id * (blocks_per_superblock - 1) + (block_id_in_superblock + !cache - 1)) + v;
		// res_upper += cache * m_block[access2];
		// res_lower += cache * m_block[access2 - 1];
		auto it2 = m_block.begin() + (access2 - 1);
		res_lower += cache * (*it2);
		++it2;
		res_upper += cache * (*it2);

        if constexpr (words_per_block > 1)
        {
            size_type const word_id{idx / values_per_word};
            uint64_t w{word_id - (word_id % words_per_block)};
            while (w < word_id)
            {
                res_upper += this->full_word_prefix_rank(this->m_v->data(), w, v);
                res_lower += this->full_word_prefix_rank(this->m_v->data(), w, v - 1);
                ++w;
            }
        }

		size_type const cache2 = idx % values_per_block != 0;
		res_upper += cache2 * this->word_prefix_rank(this->m_v->data(), idx, v);
		res_lower += cache2 * this->word_prefix_rank(this->m_v->data(), idx, v - 1);

		return res_upper - res_lower;
	}

	size_type size() const { return this->m_v->size(); }

	size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, const std::string name = "") const
	{
		structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes = 0;
		written_bytes += m_block.serialize(out, child, "prefix_block_counts");
		written_bytes += m_superblock.serialize(out, child, "prefix_superblock_counts");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	void load(std::istream& in, const int_vector<vector_width>* v = nullptr)
	{
		this->m_v = v;
		m_block.load(in);
		m_superblock.load(in);
		this->init(v);
	}

	template <typename archive_t>
	void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
	{
		ar(CEREAL_NVP(m_block));
		ar(CEREAL_NVP(m_superblock));
	}

	template <typename archive_t>
	void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
	{
		ar(CEREAL_NVP(m_block));
		ar(CEREAL_NVP(m_superblock));
	}

	void set_vector(const int_vector<vector_width>* v = nullptr)
	{
		this->m_v = v;
		this->init(v);
	}
};

} // end namespace epr
} // end namespace sdsl

#endif // end file
