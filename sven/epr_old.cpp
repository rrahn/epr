#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <numeric>

#include <iostream>
#include <seqan/align.h>
#include <typeinfo>
#include <set>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

// #include "index_fm_rank_dictionary_levels.h"


// using namespace seqan;
using namespace std;

template <typename T>
void printv(std::vector<T> & a){
    for(int i = 0; i < a.size(); ++i){
        std::cout << (int)a[i] <<",";
    }
    std::cout << "\n";
}

namespace myFunctions
{

template<typename T>
const T& mmin(const T& a, const T& b);


template <typename TSize, typename TValue, size_t TAlphabetSize, size_t TValueBitsSize, size_t levels, size_t wordsPerBlock>
struct RankDictionary
{
    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------

    static constexpr unsigned _BITS_PER_VALUE   = TValueBitsSize; //2

//     static const unsigned VALUE = sizeof(TValue) * 8;
    static constexpr unsigned _RankDictionaryWordSize = 64;
    static constexpr unsigned _BITS_PER_BLOCK   = (wordsPerBlock == 0 ? 1/*BitsPerValue<typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type>::VALUE*/ : _RankDictionaryWordSize * wordsPerBlock);

    static constexpr unsigned _BITS_PER_WORD    = mmin(_RankDictionaryWordSize, _BITS_PER_BLOCK);//Min<_RankDictionaryWordSize, _BITS_PER_BLOCK>; //64

    static constexpr unsigned _VALUES_PER_WORD  = _BITS_PER_WORD / _BITS_PER_VALUE;
    static constexpr unsigned _WORDS_PER_BLOCK  = _BITS_PER_BLOCK / _BITS_PER_WORD;
    static constexpr unsigned _VALUES_PER_BLOCK = _VALUES_PER_WORD * _WORDS_PER_BLOCK;

//     static constexpr uint64_t _VALUES_PER_SUPERBLOCK = Max<1, (((1ull << (BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 2>::Type>::VALUE/2)) - 1) / _VALUES_PER_BLOCK) * _VALUES_PER_BLOCK>::VALUE; // 2^x - 1 values
//     static constexpr uint64_t _VALUES_PER_ULTRABLOCK = Max<1, (((1ull << (BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 3>::Type>::VALUE/2)) - 1) / Max<_VALUES_PER_SUPERBLOCK, 1>::VALUE) * _VALUES_PER_SUPERBLOCK>::VALUE; // MAX: workaround for intel & clang: division is not constexpr since divisor could be 0
/*
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type   Type;

    static_assert(
        BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type>::VALUE >= LogN<_VALUES_PER_BLOCK + 1, 2>::VALUE,
        "The data type of the lowest level has to be larger (or the number of words per block smaller).");
*/

    //_RankDictionaryWordSize::TYPE
    typedef uint64_t TWordType;

    static TWordType _SELECT_BITMASK; // select only every 2nd element

    // filter by character
    static TWordType _CHAR_BITMASKS[TAlphabetSize];

    // truncate the last values in a word that shall not be counted
    static TWordType _TRUNC_BITMASKS[_VALUES_PER_WORD];

/*
    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreUltraBlocks>::Type  ultrablocks;
    typename Fibre<RankDictionary, FibreSuperBlocks>::Type  superblocks;
    typename Fibre<RankDictionary, FibreRanks>::Type        blocks;
    typename Size<RankDictionary>::Type                     _length;
    // TODO(esiragusa): open/save _length or remove it.

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    void _populateBitmasks()
    {
        // TODO(cpockrandt): make all bitmasks a const-expr
        constexpr unsigned maxValue = (1 << _BITS_PER_VALUE) - 1;
        constexpr unsigned padding = _BITS_PER_WORD % _BITS_PER_VALUE;
        // max() workaround for windows C4293 warning
        TWordType paddingWord = (padding > 0) ? (1ull << (_BITS_PER_WORD - std::max(padding, 1u))) : 0ull;
        if (isLevelsPrefixRD<Levels<TSpec, TConfig> >::VALUE && !std::is_same<TValue, bool>::value)
        {
            if (padding > 0)
                _SELECT_BITMASK = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue, 0u);
            else
                _SELECT_BITMASK = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, 0u, maxValue);
            for (unsigned i = 0; i < RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE; ++i)
            {
                if (padding > 0)
                    _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, i, 1)
                                        | paddingWord;
                else
                    _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, 1, i);
            }
            _TRUNC_BITMASKS[0] = 0u;
            for (unsigned i = 0;
                 i < RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >::VALUE - 1;
                 ++i)
            {
                if (padding > 0)
                    _TRUNC_BITMASKS[i+1] = _bitmask<TWordType>(_BITS_PER_WORD,2*i + 1, _BITS_PER_VALUE, 0, 1)
                                           | paddingWord;
                else
                    _TRUNC_BITMASKS[i+1] = _bitmask<TWordType>(_BITS_PER_WORD, 2*i + 1, _BITS_PER_VALUE, 1, 0);
            }
        }
        else
        {
            for (unsigned i = 0; i < RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE; ++i)
                _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD,
                                                        _VALUES_PER_WORD,
                                                        _BITS_PER_VALUE,
                                                        maxValue-i,
                                                        maxValue-i);
            for (unsigned i = 0;
                 i < RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >::VALUE;
                 ++i)
                _TRUNC_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, i+1, _BITS_PER_VALUE, 1, 1);
        }
    }

    RankDictionary() :
        _length(0)
    {
        _populateBitmasks();
    }

    template <typename TText>
    RankDictionary(TText const & text) :
        _length(0)
    {
        _populateBitmasks();

        createRankDictionary(*this, text);
    }*/

};


// #define TRD RankDictionary<TValue, Levels<TSpec, TConfig> >
#define TRD RankDictionary<TSize, TValue, TAlphabetSize, TValueBitsSize, levels, wordsPerBlock>
template <typename TSize, typename TValue, size_t TAlphabetSize, size_t TValueBitsSize, size_t levels, size_t wordsPerBlock>
typename TRD::TWordType TRD::_SELECT_BITMASK;

template <typename TSize, typename TValue, size_t TAlphabetSize, size_t TValueBitsSize, size_t levels, size_t wordsPerBlock>
typename TRD::TWordType TRD::_CHAR_BITMASKS[TAlphabetSize];

template <typename TSize, typename TValue, size_t TAlphabetSize, size_t TValueBitsSize, size_t levels, size_t wordsPerBlock>
typename TRD::TWordType TRD::_TRUNC_BITMASKS[_VALUES_PER_WORD];
#undef TRD


template <typename T>
int lengthSum(T & a){
    unsigned lengthSum{0};
    for(unsigned i = 0; i < a.size(); ++i)
        lengthSum += length(a[i]);
    return(lengthSum);
}


}

template <typename TBwt, typename TSentinels, typename TText, typename TSA>
inline void
createBwt(TBwt & bwt, TSentinels & sentinels, TText const & text, TSA const & sa)
{
    typedef typename seqan::Value<TSA>::Type                       TSAValue;
    typedef typename seqan::Size<TSA>::Type                        TSize;

    // Resize the RankDictionary.
    TSize seqNum = text.size();
    TSize totalLen = seqan::length(sa); //TODO why is length sum to short?? //myFunctions::lengthSum(text);

    std::cout << "seqNum " << seqNum << "\ttotalLen " << totalLen << "\tSAlength " << length(sa) << "\n";

    int sentinelSubstitute = 0;
//     std::vector<bool> sentinels(seqNum + totalLen, false);



    resize(sentinels, seqNum + totalLen, seqan::Exact());
    resize(bwt, seqNum + totalLen, seqan::Exact());

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 0; i < seqNum; ++i)
    {
        //TODO change to size
        if (length(text[seqNum - (i + 1)]) > 0)
        {
            bwt[i] =  back(text[seqNum - (i + 1)]);
            sentinels[i] = false;
        }
    }

    // Compute the rest of the bwt.
    for (TSize i = 0; i < length(sa); ++i)
    {
        TSAValue pos = sa[i];
//         posLocalize(pos, sa[i], stringSetLimits(text));

        if (getSeqOffset(pos) != 0)
        {
//             std::cout << i << " " << text[getSeqNo(pos)][getSeqOffset(pos) - 1] << " "<< (ordValue(text[getSeqNo(pos)][getSeqOffset(pos) - 1]) + 1) << "\n";
            bwt[i + seqNum] = ordValue(text[getSeqNo(pos)][getSeqOffset(pos) - 1]) + 1;
            sentinels[i + seqNum] = false;
        }
        else
        {
            bwt[i + seqNum] = sentinelSubstitute;
            sentinels[i + seqNum] = true;
        }
    }
}


template <typename TChar, typename TConfig>
void generateText(seqan::String<TChar, TConfig> & text, unsigned const length)
{
    unsigned alphabetSize = seqan::ValueSize<TChar>::VALUE;
    unsigned minChar = seqan::MinValue<TChar>::VALUE;

    seqan::resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = std::rand() % alphabetSize - minChar;
    }
}

template <typename TPrefixSums, typename TText, typename TValue>
inline void prefixSums(TPrefixSums & sums, TText const & text, TValue alphabetSize)
{
    typedef typename seqan::Concatenator<TText const>::Type        TConcat;
    typedef typename seqan::Iterator<TConcat, seqan::Standard>::Type      TIter;

    seqan::resize(sums, alphabetSize + 1, 0, seqan::Exact());

    // Compute symbol frequencies.
    TIter itEnd = seqan::end(seqan::concat(text), seqan::Standard());
    for (TIter it = seqan::begin(seqan::concat(text), seqan::Standard()); it != itEnd; seqan::goNext(it))
        sums[seqan::ordValue(*it) + 1]++;

    // Cumulate symbol frequencies.
    for(int i = 1; i < sums.size(); ++i)
        sums[i] += sums[i - 1];
}



int main(int argc, char const * argv[])
{
    seqan::ArgumentParser parser("Test");

    seqan::addOption(parser, seqan::ArgParseOption("l", "length", "length of the read", seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption("e", "errors", "Number of allowed Errors", seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption("m", "mutations", "Number of mutations done on the read", seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption("i", "iterations", "Number of iterations testing search", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("HD", "hammingDistance", ""));

    seqan::addOption(parser, seqan::ArgParseOption("em", "editMuations", "Allow Insertions and Deletions as mutations"));

    seqan::addOption(parser, seqan::ArgParseOption("mv", "myversion", "Use my version used in yara"));

//     addOption(parser, ArgParseOption("iv", "itv", "Use my version with in text verification"));

    seqan::addOption(parser, seqan::ArgParseOption("v", "verbose", ""));

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;


    int textlength = 12;

    typedef seqan::Dna     TAlphabet;
    typedef seqan::DnaString TText;
    typedef seqan::StringSet<TText> TSeqan2Text;
    typedef std::vector<TText> TSeqan3Text;
    typedef seqan::FastFMIndexConfig<void, uint32_t> TMyFastConfig;
    typedef seqan::Index<TSeqan2Text, seqan::FMIndex<void, TMyFastConfig> > TIndex;
    typedef seqan::Iter<TIndex, seqan::VSTree<seqan::TopDown<> > > TIter;

    TSeqan2Text text;
    TSeqan3Text cText;
    std::vector<int> sums;
    std::vector<uint8_t> bwt;
    std::vector<bool> sentinels;

    TText rand;
    int i = 0;
    std::cout << "Text:\n";
    while(i < 3){
        generateText(rand, textlength);
        std::cout << rand << "\n";
        seqan::appendValue(text, rand);
        cText.push_back(rand);
        ++i;
    }
    std::cout << "\n";

    TIndex index(text);
    TIter it(index);





    unsigned alphSize = seqan::ValueSize<TAlphabet>::VALUE;
    prefixSums(sums, text, alphSize);

    unsigned sentinelsCount = cText.size();
    for (unsigned i = 0; i < seqan::length(sums); ++i)
        sums[i] += sentinelsCount;


    createBwt(bwt, sentinels, cText, index.sa);

    std::cout << "bwt\n";
    printv(bwt);
    std::cout << "sen\n";
    printv(sentinels);

// #include <cmath>
//     std::log2(0.125)

    //TSize, TValue, TValueSize, TLevels, WordsPerBlock
    myFunctions::RankDictionary<uint32_t, uint8_t, 4, 2, 2, 1> bwts;
//     = updateRanks(bwt);

/*
   // Update all ranks.
   updateRanks(lf.bwt);
   // Update the auxiliary RankDictionary of sentinel positions.
   updateRanks(lf.sentinels);*/

/*
    uint32_t indexSize = seqan::length(index.sa);
    std::cout << "IndexSize: " << indexSize << "\n";

    std::cout << "Print sums: \n";
    printv(sums);

    for(int i = 0; i < indexSize; ++i){
        std::cout << i << ":\t" << index.sa[i] << "\n";
    }*/



    return 0;
}



