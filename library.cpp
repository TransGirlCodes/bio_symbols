#include "library.h"

#include <iostream>

template <typename NucleicAcid>
NucleicAcid operator&(NucleicAcid a, NucleicAcid b){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<DNA>::type;
    return static_cast<NucleicAcid>(static_cast<UT>(a) & static_cast<UT>(b));
}

template <typename NucleicAcid>
NucleicAcid operator|(NucleicAcid a, NucleicAcid b){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<NucleicAcid>(static_cast<UT>(a) | static_cast<UT>(b));
}

template <typename NucleicAcid>
NucleicAcid operator~(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<DNA>(~static_cast<UT>(nt) & 0b1111);
}

template<typename NucleicAcid>
bool isGC(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::W) == NucleicAcid::Gap;
}

template<typename NucleicAcid>
bool ispurine(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::Y) == NucleicAcid::Gap;
}

template<typename NucleicAcid>
bool ispyrimidine(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    return nt != NucleicAcid::Gap and (nt & NucleicAcid::R) == NucleicAcid::Gap;
}

bool isambiguous(DNA nt){
    // need popcount
}

bool iscertain(DNA nt){
    // need popcount
}

bool isgap(DNA nt){
    // need popcount
}

template <typename NucleicAcid>
NucleicAcid complement(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    UT bits = static_cast<UT>(nt);
    return static_cast<NucleicAcid>((bits & 0b0001) << 3 | (bits & 0b1000) >> 3 |
                                    (bits & 0b0010) << 1 | (bits & 0b0100) >> 1);
}

bool isvalid(int x){
    return 0 <= x and x < 16;
}

template <typename NucleicAcid>
bool isvalid(NucleicAcid nt){
    static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids");
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return static_cast<UT>(nt) <= 0b1111;
}
/*
function Base.isvalid{T<:NucleicAcid}(::Type{T}, x::Integer)
    return 0 ≤ x < 16
end

function Base.isvalid(nt::NucleicAcid)
    return reinterpret(UInt8, nt) ≤ 0b1111
end

@inline function iscompatible{T<:NucleicAcid}(x::T, y::T)
    return compatbits(x) & compatbits(y) != 0
end
*/

// Test if `x` and `y` are compatible with each other (i.e. `x` and `y` can be the same symbol).
template <typename NucleicAcid>
bool iscompatible(NucleicAcid x, NucleicAcid y){
    using UT = typename std::underlying_type<NucleicAcid>::type;
    return (static_cast<UT>(x) & static_cast<UT>(y)) != 0;
}
