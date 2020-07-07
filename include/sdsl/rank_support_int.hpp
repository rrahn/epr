// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int.hpp
    \brief rank_support_int.hpp contains classes that support a sdsl::int_vector with constant time rank information.
           Rank is defined as the number of occurrences of a value up to a given position.
	\author Christopher Pockrandt
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT
#define INCLUDED_SDSL_RANK_SUPPORT_INT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::int_vector with the rank method.
 */

#include "int_vector.hpp"
#include "uint128_t.hpp"

#include <bitset>

// TODO: benchmark the use of compiler hints for branch prediction
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

//! Namespace for the succinct data structure library.
namespace sdsl {
namespace epr {

//! The base class of classes supporting rank_queries for a sdsl::int_vector in constant time.
/*!
*/

constexpr size_t floor_log2(size_t const n)
{
    return (n == 1) ? 0 : 1 + floor_log2(n >> 1);
}

constexpr size_t ceil_log2(size_t const n)
{
    return (n == 1) ? 0 : floor_log2(n - 1) + 1;
}

template <uint8_t alphabet_size, uint8_t vector_width>
class rank_support_int {

public:
	typedef typename int_vector<vector_width>::size_type size_type;
	typedef typename int_vector<vector_width>::value_type value_type;

    static_assert(alphabet_size > 2, "Rank support is only implemented on int_vectors with an alphabet size of > 2.");

protected:

    // Constructs a bit mask with the pattern w of a given length.
    // It is concatenated until the length of the bitmask reaches max_length.
    template <typename uintX_t>
    static constexpr uintX_t bm_rec(const uintX_t w, const uint8_t length, const uint8_t max_length)
    {
        return (length >= max_length) ? w : bm_rec(w | (w << length), length << 1, max_length);
    }

protected:
	const int_vector<vector_width>* m_v; //!< Pointer to the rank supported bit_vector
    const uint64_t * data_ptr;
    static constexpr uint8_t sigma{alphabet_size};
    // static constexpr uint8_t values_per_word{64/sigma};
    static constexpr uint8_t sigma_bits{ceil_log2(alphabet_size)};
    static constexpr uint8_t bits_per_word{(64 / sigma_bits) * sigma_bits};
    // static constexpr uint8_t sigma_bits{alphabet_size};
	static constexpr uint64_t even_mask{bm_rec<uint64_t>(bits::lo_set[sigma_bits], sigma_bits * 2, 64)};
	static constexpr uint64_t carry_select_mask{bm_rec<uint64_t>(1ULL << sigma_bits, sigma_bits * 2, 64)};
	uint64_t masks[alphabet_size]; // TODO: make constexpr and remove init()

    // Count how often value v or smaller occurs in the word w.
    uint64_t set_positions_prefix(const uint64_t w, const value_type v) const
    {
        uint64_t const w_even = even_mask & w; // retrieve even positions
        uint64_t const w_odd = even_mask & (w >> sigma_bits); // retrieve odd positions
        // every bit that will be set corresponds to an element <= v
        // because the preset bit to the left in the precomputed bitmask is not eliminated by the carry bit during the subtraction
        uint64_t res = (masks[v] - w_even) & carry_select_mask;
        // since alphabet_size is > 2 and an element uses at least 2 bits, we can shift the odd positions by one to the left
        // and it is guaranteed that when adding both with OR that no bits that are set will overlap.
        res |= ((masks[v] - w_odd) & carry_select_mask) << 1;
        return res;
    }

    // Count how often value v occurs in the word w.
    // Cannot be called on v = 0. Call set_positions_prefix(w, 0) instead.
    uint64_t set_positions(const uint64_t w, const value_type v) const
    {
        assert(v > 0);
        // optimiyed version of set_positions(w, v) - set_positions(w, v - 1)
        uint64_t const w_even = even_mask & w; // retrieve even positions
        uint64_t const w_odd = even_mask & (w >> sigma_bits); // retrieve odd positions
        uint64_t res = ((masks[v] - w_even) & ~(masks[v - 1] - w_even)) & carry_select_mask;
        res |= (((masks[v] - w_odd) & ~(masks[v - 1] - w_odd)) & carry_select_mask) << 1;
        return res;
    }

    // Counts the occurrences of elements smaller or equal to v in the word starting at data up to position idx.
    uint32_t word_prefix_rank_(const uint64_t word, const size_type bit_pos, const value_type v) const
    {
        // return v;
        // std::cout << "bit_pos=" << bit_pos << '\n';
        // std::cout << "w= " << std::bitset<64>(w).to_string() << " v=" << std::bitset<64>(v).to_string() << " set_positions_prefix=" << std::bitset<64>(set_positions_prefix(w, v)).to_string() <<
        //              " bits::lo_set=" << std::bitset<64>(bits::lo_set[(bit_pos % 63)]).to_string() << '\n';
        return bits::cnt(set_positions_prefix(word, v) & bits::lo_set[(bit_pos % bits_per_word) + 1] ); //bits::lo_set[(bit_pos & 0x3F) + 1 + bit_pos/63]
    }

    // Counts the occurrences of elements smaller or equal to v in the word starting at data up to position idx.
    uint32_t word_prefix_rank(const uint64_t* data, const size_type idx, const value_type v) const
    {
        size_type const bit_pos = idx * sigma_bits;
        return word_prefix_rank_(*(data + (bit_pos / bits_per_word)), bit_pos, v);
        // uint64_t const w = ;
        // return v;
        // std::cout << "bit_pos=" << bit_pos << '\n';
        // std::cout << "w= " << std::bitset<64>(w).to_string() << " v=" << std::bitset<64>(v).to_string() << " set_positions_prefix=" << std::bitset<64>(set_positions_prefix(w, v)).to_string() <<
        //              " bits::lo_set=" << std::bitset<64>(bits::lo_set[(bit_pos % 63)]).to_string() << '\n';
        // return bits::cnt(set_positions_prefix(w, v) & bits::lo_set[(bit_pos % bits_per_word) + 1] ); //bits::lo_set[(bit_pos & 0x3F) + 1 + bit_pos/63]
    }

    // Counts the occurrences of elements smaller or equal to v in the word starting at data up to position idx.
    // Cannot be called on v = 0. Call word_prefix_rank(data, idx, 0) instead.
    uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v) const
    {
        size_type const bit_pos = idx * sigma_bits;
        uint64_t const w = *(data + (bit_pos >> 6));
        return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
    }

    // Counts the occurrences of v in the word starting at data up to position idx.
    uint32_t full_word_prefix_rank(const uint64_t* data, const size_type word_pos, const value_type v) const
    {
        uint64_t const w = *(data + word_pos);
        return bits::cnt(set_positions_prefix(w, v));
    }

    // Counts the occurrences of v in the word starting at data up to position idx.
    // Cannot be called on v = 0. Call full_word_prefix_rank(data, word_pos, 0) instead.
    uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v) const
    {
        uint64_t const w = *(data + word_pos);
        return bits::cnt(set_positions(w, v));
    }

    // Constructs a bitmask for each value of the alphabet.
    void init(const int_vector<vector_width>* v)
    {
        if (v != nullptr) {
            assert(sigma_bits == v->width()); // TODO Not valid because EPR uses effective sigma, the text uses sigma
            m_v = v;
            data_ptr = m_v->data();

            for (value_type v = 0; v < alphabet_size; ++v)
            {
                masks[v] = v;
                for (uint8_t i = sigma_bits * 2; i < 64; i <<= 1)
                    masks[v] |= masks[v] << i;
            }

            uint64_t tmp_carry = masks[1];
            for (value_type v = 0; v < alphabet_size; ++v)
                masks[v] |= tmp_carry << sigma_bits;
        }
    }

public:
	//! Constructor
	/*! \param v The supported int_vector.
         */
	rank_support_int(const int_vector<vector_width>* v = nullptr);
	//! Copy constructor
	rank_support_int(const rank_support_int&) = default;
	rank_support_int(rank_support_int&&)	  = default;
	rank_support_int& operator=(const rank_support_int&) = default;
	rank_support_int& operator=(rank_support_int&&) = default;
	//! Destructor
	virtual ~rank_support_int() {}

	//! Answers rank queries for the supported int_vector.
	/*!	\param i Argument for the length of the prefix v[0..i-1].
            \param v Argument which value to count.
        	\returns Number of occurrences of v in the prefix [0..i-1] of the supported int_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
	virtual size_type rank(const size_type i, const value_type v) const = 0;
	//! Alias for rank(idx, v)
	virtual size_type operator()(const size_type idx, const value_type v) const = 0;
	//! Answers rank queries for the supported int_vector.
	/*!	\param i Argument for the length of the prefix v[0..i-1].
            \param v Argument which value (including smaller values) to count.
        	\returns Number of occurrences of elements smaller or equal to v in the prefix [0..i-1] of the supported int_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
	virtual size_type prefix_rank(const size_type i, const value_type v) const = 0;
	//! Serializes rank_support_int.
	/*! \param out Out-Stream to serialize the data to.
        */
	virtual size_type
	serialize(std::ostream& out, structure_tree_node* v, const std::string name) const = 0;
	//! Loads the rank_support_int.
	/*! \param in In-Stream to load the rank_support_int data from.
            \param v The supported int_vector.
         */
	virtual void load(std::istream& in, const int_vector<vector_width>* v = nullptr) = 0;
	//! Sets the supported int_vector to the given pointer.
	/*! \param v The new int_vector to support.
         *  \note Method init has to be called before the next call of rank or prefix_rank.
         *  \sa init, rank, prefix_rank
         */
	virtual void set_vector(const int_vector<vector_width>* v = nullptr) = 0;
};

template <uint8_t alphabet_type, uint8_t vector_width>
inline rank_support_int<alphabet_type,vector_width>::rank_support_int(const int_vector<vector_width>* v)
{
    init(v);
}

} // end namespace epr
} // end namespace sdsl

#include "rank_support_int_v.hpp"
#include "rank_support_int_scan.hpp"

#endif // end file
