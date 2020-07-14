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

#include <sdsl/cereal.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl {
namespace epr {


template <typename value_t, size_t bits_per_value>
class bit_compressed_word
{
private:
    static constexpr uint64_t max_size = (sizeof(uint64_t) << 3) / bits_per_value;
    static constexpr uint64_t bit_mask = bits::lo_set[bits_per_value];
    static constexpr uint64_t mask = ((1ull << (max_size * bits_per_value - 1)) - 1ull) << 1 | 1ull;

    uint64_t word{};

public:

    bit_compressed_word() = default;

    template <typename it_t>
    constexpr bit_compressed_word(it_t it, it_t end) noexcept
    {
        assign(it, end);
    }

    constexpr value_t operator[](size_t const index) const noexcept
    {
        assert(index < max_size);
        uint64_t offset = index * bits_per_value;
        return value_t{(word >> offset) & bit_mask};
    }

    template <typename it_t>
    constexpr void assign(it_t it, it_t end) noexcept
    {
        assert(std::distance(it, end) <= max_size);

        for (size_t index = 0; it != end; ++it, ++index)
        {
            uint64_t offset = index * bits_per_value;
            word = (word & ~(bit_mask << offset)) | uint64_t{*it} << offset;
        }
    }

    // template <typename TPos, typename TValue2>
    // inline TValue2
    // assignValue(TPos k, TValue2 const source)
    // {
    //     SEQAN_ASSERT_GEQ(static_cast<int64_t>(k), 0);
    //     SEQAN_ASSERT_LT(static_cast<int64_t>(k), static_cast<int64_t>(SIZE));

    //     unsigned shift = ((SIZE - 1 - k) * BitsPerValue<TValue>::VALUE);
    //     i = (i & ~(BIT_MASK << shift)) | (TBitVector)ordValue(source) << shift;
    //     return source;
    // }

    constexpr std::add_pointer_t<uint64_t>  data() noexcept
    {
        return &word;
    }

    constexpr std::add_pointer_t<uint64_t const> data() const noexcept
    {
        return &word;
    }

    constexpr operator uint64_t() const noexcept
    {
        return word;
    }
};

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
class rank_support_int_v final : public rank_support_int<alphabet_size,vector_width> {
public:
	typedef typename rank_support_int<alphabet_size,vector_width>::size_type size_type;
	typedef typename rank_support_int<alphabet_size,vector_width>::value_type value_type;

	using rank_support_int<alphabet_size,vector_width>::sigma;
	using rank_support_int<alphabet_size,vector_width>::sigma_bits;

public:
	static constexpr uint64_t values_per_word{64ULL / sigma_bits};
	static constexpr uint32_t values_per_block{words_per_block * values_per_word};
    static constexpr uint64_t words_per_superblock{words_per_block * blocks_per_superblock};
    static constexpr uint64_t values_per_superblock{blocks_per_superblock * values_per_block};

private:

	struct superblock_node;

    std::vector<superblock_node> nodes{};

public:
	explicit rank_support_int_v(const int_vector<vector_width>* v = nullptr) : rank_support_int<alphabet_size,vector_width>(v)
	{
	    static_assert(blocks_per_superblock > 1, "There must be at least two blocks per superblock!");
		constexpr size_t max_letter{sigma - 1};

		if (v == nullptr)
		{
			return;
		}
		else if (v->empty())
		{
			// m_block.resize(max_letter, 0);
			// m_superblock.resize(max_letter, 0);
			return;
		}

		constexpr uint64_t new_width{ceil_log2(values_per_superblock)};

		// NOTE: number of elements is artificially increased by one because rank can be called on m_v[size()]
		uint64_t const word_count = ((this->m_v->size() - 1 + 1) / values_per_word) + 1; // equivalent to ceil(m_v->size() / values_per_word)
		uint64_t const block_count = ((word_count - 1) / words_per_block) + 1; // equivalent to ceil(word_count / words_per_block)

		// for each superblock we only need `blocks_per_superblock-1` instead of `blocks_per_superblock` blocks.
		// for the last superblock we can subtract the last unused blocks.
        size_type const blocks_needed = (((block_count - 1) / blocks_per_superblock) + 1) * (blocks_per_superblock - 1)
									  - ((blocks_per_superblock - (block_count % blocks_per_superblock)) % blocks_per_superblock);
		size_type const block_size = (blocks_needed + 1) * max_letter;
		size_type const superblock_count = (word_count + words_per_superblock - 1) / words_per_superblock;
		size_type const superblock_size = superblock_count * max_letter; // equivalent to ceil(word_count / words_per_superblock) * max_letter

		uint64_t const * data = this->m_v->data();
		std::array<uint64_t, max_letter> buf_blocks{};
		std::array<uint64_t, max_letter> buf_superblocks{};

        // Allocate a node for every superblock.
        nodes.resize(superblock_count);

        // Initialise local buffers to keep track of the cumulative sum.
		buf_blocks.fill(0);
		buf_superblocks.fill(0);

		// Precompute blocks and superblocks
		auto text_slice_it = v->begin();
        uint64_t word_id = 0;  // We basically iterate over all words of the underlying text.
        // auto superblock_it = m_superblock.begin();
        size_t superblock_id = 0;
        size_t block_id = 0;
        // auto block_it = m_block.begin();
        for (auto node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
            // First initialise the superblock text slice.
            for (auto & compressed_word : node_it->superblock_text)
            {
                auto text_slice_end = std::next(text_slice_it,
                                                std::min<size_t>(std::distance(text_slice_it, v->end()),
                                                                 values_per_word));
                compressed_word.assign(text_slice_it, text_slice_end);
                text_slice_it = text_slice_end;
            }

            // Second initialise the superblock counts: buf_blocks corresponds to
            auto superblock_it = node_it->superblocks.begin();  // Store the beign of the super block in the node.
            for (size_t letter_rank = 0; letter_rank < max_letter; ++letter_rank, ++superblock_it, ++superblock_id)
            {
                buf_superblocks[letter_rank] += buf_blocks[letter_rank]; // Update sum with previous superblock
                *superblock_it = buf_superblocks[letter_rank]; // Store the counts.
                buf_blocks[letter_rank] = 0; // Reset the block counts for the next superblock.
            }

            // Third initialise the block counts:
            // The stored block counts represent the cumulative sum of the previous blocks in the super block.
            // The first block of the superblock is not stored explicitly since it has no predecessor.
            // A block stores the counts for the letters consecutive in memory from [0..max_letter] and starts then the
            // next block at offset i*max_letter, where i is the current block id.
            // TODO: Make the implementation safe for multiple words per block
            for (auto block_it = node_it->blocks.begin(); word_id < word_count && block_it != node_it->blocks.end(); ++word_id)
            {
                // Get the prefix ranks for the current word for each letter and store them in the respective block
                for (size_t letter_rank = 0; letter_rank < max_letter; ++letter_rank, ++block_it, ++block_id)
                {
                    buf_blocks[letter_rank] += this->full_word_prefix_rank(data, word_id, letter_rank);
                    *block_it = buf_blocks[letter_rank];
                }
            }

            // Count the last block which is not stored explicitly.
            if (word_id < word_count)
            {
                for (uint64_t letter = 0; letter < max_letter; ++letter)
                    buf_blocks[letter] += this->full_word_prefix_rank(data, word_id, letter);

                ++word_id;
            }
        }

        // Compare the blocks: They should store the same values.

        // auto cmp_it = m_block.begin();
        // size_t index = 0;
        // size_t node_id = 0;
        // for (auto node_it = nodes.begin(); node_it != nodes.end() && cmp_it != m_block.end(); ++node_it)
        // {
        //     auto && node = *node_it;
        //     size_t block_id = 0;
        //     auto it = node.blocks.begin();
        //     for (size_t i = 0; i < ((blocks_per_superblock - 1) * max_letter) && cmp_it != m_block.end(); ++i, ++it)
        //     {
        //         if (*cmp_it != *it)
        //             std::cout << "FAIL: m_block[" << std::setw(4) << index << "] != node["
        //                       << std::setw(3) << node_id << "]blocks[" << std::setw(3) << block_id << "]: "
        //                       << std::setw(2) << *cmp_it << " != " << std::setw(2) << *it << "\n";

        //         ++index;
        //         ++block_id;
        //         ++cmp_it;
        //     }
        //     ++node_id;
        // }
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

        size_type const node_id = idx / values_per_superblock;
        size_type const block_id_in_superblock = (idx / values_per_block) % blocks_per_superblock;

        auto const & node = nodes[node_id];

		// retrieve superblock value
        size_type res = node.superblocks[v];

		// retrieve block value
        // The first block stores the counts for the actual second block.
        // Hence, if we are in the first block of the superblock we get the first block value but multiply it with
        // 0 (cache * node.block[pppp]). If the block id is greater than 0 we remove one from the letter
		size_type const is_first_block = block_id_in_superblock == 0;
		size_type const pppp = (max_letter * (block_id_in_superblock + is_first_block - 1)) + v;
		res += !is_first_block * (node.blocks[pppp]);

        // TODO: Enable me!
		// compute in-block queries for all words before the in-block queries
		// this only applies when multiple words are in one block
        // if constexpr (words_per_block > 1)
        // {
        //     size_type const word_id{idx / values_per_word};
        //     uint64_t w{word_id - (word_id % words_per_block)};
        //     while (w < word_id)
        //     {
        //         res += this->full_word_prefix_rank(this->m_v->data(), w, v);
        //         ++w;
        //     }
		// 	// std::cout << "res3=" << res << '\n';
        // }

		// compute in-block query
		size_type cache2 = idx % values_per_block != 0;
		// size_type const cache2 = 1;
        // First, get the local bit position within the data of the super block.
        size_type const bit_pos = (idx % values_per_superblock) * sigma_bits;
        // Second, compute the word that contains this value.
        uint64_t const word = node.superblock_text[bit_pos / this->bits_per_word];
        // Third, compute the in-block rank given the current word.
        res += cache2 * this->word_prefix_rank_(word, bit_pos, v);

		return res;
	}

	size_type prefix_rank2(const size_type idx, const value_type v) const
	{
		assert(this->m_v != nullptr);
		assert(idx <= this->m_v->size());
        assert(v < this->sigma); // v cannot have rank sigma - 1
        assert(v > 0); // v cannot have rank 0

		constexpr uint8_t const max_letter{this->sigma - 1};

        size_type const node_id = idx / values_per_superblock;
        size_type const block_id_in_superblock_ = (idx / values_per_block) % blocks_per_superblock;
        auto && node = nodes[node_id];

        // First access is the superblock_id offset with the max_letter plus the current letter.
        size_t res_lower_ = node.superblocks[v - 1];
        size_t res_upper_ = node.superblocks[v];

		size_type const is_first_block = block_id_in_superblock_ == 0;

        // Now we do the same with the blocks. We know that v is at least 1.
        size_type pppp = (max_letter * (block_id_in_superblock_ + is_first_block - 1)) + (v - 1);
        res_lower_ += !is_first_block * node.blocks[pppp];
        res_upper_ += !is_first_block * node.blocks[++pppp];

        // TODO: Enable me!
        // if constexpr (words_per_block > 1)
        // {
        //     size_type const word_id{idx / values_per_word};
        //     uint64_t w{word_id - (word_id % words_per_block)};
        //     while (w < word_id)
        //     {
        //         res_upper += this->full_word_prefix_rank(this->m_v->data(), w, v);
        //         res_lower += this->full_word_prefix_rank(this->m_v->data(), w, v - 1);
        //         ++w;
        //     }
        // }

        size_type cache2_ = idx % values_per_block != 0;
		// size_type const cache2 = 1;
        // First, get the local bit position within the data of the super block.
        size_type const bit_pos_ = (idx % values_per_superblock) * sigma_bits;
        // Second, compute the word that contains this value.
        uint64_t const word_ = node.superblock_text[bit_pos_ / this->bits_per_word];
        // Third, compute the in-block rank given the current word.
        res_lower_ += cache2_ * this->word_prefix_rank_(word_, bit_pos_, v - 1);
        res_upper_ += cache2_ * this->word_prefix_rank_(word_, bit_pos_, v);

		return res_upper_ - res_lower_;
	}

	size_type size() const { return this->m_v->size(); }

	size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, const std::string name = "") const
	{
		structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes = 0;
        // for (auto & node : nodes)
        //     written_bytes += node.serialize(out, child, "nodes");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	void load(std::istream& in, const int_vector<vector_width>* v = nullptr)
	{
		this->m_v = v;
        // while (!in.eof())
        // {
        //     superblock_node node;
        //     node.load(in);
        //     nodes.push_back(std::move(node));
        // }
		this->init(v);
	}

	template <typename archive_t>
	void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
	{
		ar(CEREAL_NVP(nodes));
	}

	template <typename archive_t>
	void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
	{
		ar(CEREAL_NVP(nodes));
	}

	void set_vector(const int_vector<vector_width>* v = nullptr)
	{
		this->m_v = v;
		this->init(v);
	}
};

template <uint8_t alphabet_size, uint8_t vector_width, uint8_t words_per_block, uint8_t blocks_per_superblock>
struct rank_support_int_v<alphabet_size, vector_width, words_per_block, blocks_per_superblock>::superblock_node
{
    // typename int_vector<64>::const_iterator superblock_it;
    // typename int_vector<0>::const_iterator block_it;
    std::array<uint32_t, (alphabet_size - 1)> superblocks;
    std::array<uint32_t, (blocks_per_superblock - 1) * (alphabet_size - 1)> blocks;
    std::array<bit_compressed_word<uint8_t, sigma_bits>, words_per_superblock> superblock_text;
    // std::array<compressed_block_text<values_per_word>, words_per_superblock> superblock_text;

    // size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, const std::string name = "") const
	// {
	// 	structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
	// 	size_type written_bytes = 0;
	// 	written_bytes += blocks.serialize(out, child, "prefix_block_counts");
	// 	written_bytes += superblocks.serialize(out, child, "prefix_superblock_counts");
	// 	written_bytes += text_slice.serialize(out, child, "superblock_text_slice");
	// 	structure_tree::add_size(child, written_bytes);
	// 	return written_bytes;
	// }

	// void load(std::istream& in)
	// {
	// 	blocks.load(in);
	// 	superblocks.load(in);
	// 	text_slice.load(in);
	// }

	// template <typename archive_t>
	// void CEREAL_SAVE_FUNCTION_NAME(archive_t & ar) const
	// {
	// 	ar(CEREAL_NVP(blocks));
	// 	ar(CEREAL_NVP(superblocks));
	// 	ar(CEREAL_NVP(text_slice));
	// }

	// template <typename archive_t>
	// void CEREAL_LOAD_FUNCTION_NAME(archive_t & ar)
	// {
	// 	ar(CEREAL_NVP(blocks));
	// 	ar(CEREAL_NVP(superblocks));
	// 	ar(CEREAL_NVP(text_slice));
	// }
};

} // end namespace epr
} // end namespace sdsl

#endif // end file
