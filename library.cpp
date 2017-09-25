#include "library.h"
#include <type_traits>

template <typename NucleicAcid>
int count_ones(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    auto b = static_cast<UT>(nt);
#if __has_builtin(__builtin_popcount)
    return __builtin_popcount(b);
#else
    b = b - ((b >> 1) & 0x55);
    b = (b & 0x33) + ((b >> 2) & 0x33);
    return (((b + (b >> 4)) & 0x0F) * 0x01);
#endif
}

template <typename NucleicAcid>
NucleicAcid operator&(NucleicAcid a, NucleicAcid b){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    using UT = typename std::underlying_type<DNA>::type;
    return static_cast<NucleicAcid>(static_cast<UT>(a) & static_cast<UT>(b));
}

template <typename NucleicAcid>
NucleicAcid operator|(NucleicAcid a, NucleicAcid b){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<NucleicAcid>(static_cast<UT>(a) | static_cast<UT>(b));
}

template <typename NucleicAcid>
NucleicAcid operator~(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<DNA>(~static_cast<UT>(nt) & 0b1111);
}

template<typename NucleicAcid>
bool is_GC(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::W) == NucleicAcid::Gap;
}

template<typename NucleicAcid>
bool is_purine(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::Y) == NucleicAcid::Gap;
}

template<typename NucleicAcid>
bool is_pyrimidine(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::R) == NucleicAcid::Gap;
}

template <typename NucleicAcid>
bool is_ambiguous(NucleicAcid nt){
    return count_ones(nt) > 1;
}

bool is_certain(DNA nt){
    return count_ones(nt) == 1;
}

template <typename NucleicAcid>
bool is_gap(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    return nt == NucleicAcid::Gap;
}

template <typename NucleicAcid>
NucleicAcid complement(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    auto bits = static_cast<UT>(nt);
    return static_cast<NucleicAcid>((bits & 0b0001) << 3 | (bits & 0b1000) >> 3 |
                                    (bits & 0b0010) << 1 | (bits & 0b0100) >> 1);
}

bool is_valid(int x){
    return 0 <= x and x < 16;
}

template <typename NucleicAcid>
bool is_valid(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<UT>(nt) <= 0b1111;
}

// Test if `x` and `y` are compatible with each other (i.e. `x` and `y` can be the same symbol).
template <typename NucleicAcid>
bool is_compatible(NucleicAcid x, NucleicAcid y){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return (static_cast<UT>(x) & static_cast<UT>(y)) != 0;
}