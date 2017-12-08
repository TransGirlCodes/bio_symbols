#ifndef BIOSYMBOLS_H
#define BIOSYMBOLS_H

#include <type_traits>

namespace biosymbols {
    template<typename T>
    struct is_nucleic_acid{
        static const bool value = false;
    };

// The full gamut of symbols for the DNA nucleotide alphabet.
    enum class DNA : unsigned char { Gap, A, C, M, G, R, S, V,
				     T, W, Y, H, K, D, B, N, Invalid };

    template<>
    struct is_nucleic_acid<DNA>{
        static const bool value = true;
    };

// The full gamut of symbols for the RNA alphabet.
    enum class RNA : unsigned char { Gap, A, C, M, G, R, S, V,
				     U, W, Y, H, K, D, B, N, Invalid };

    template<>
    struct is_nucleic_acid<RNA>{
        static const bool value = true;
    };

// The full gamut of symbols for the amino acid alphabet.
    enum class AA : unsigned char { A, R, N, D, C, Q, E, G, H, I,
	                            L, K, M, F, P, S, T, W, Y, V,
				    O, U, B, J, Z, X, Term, Gap };

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
    NucleicAcid not_version(NucleicAcid nt){
        static_assert(is_nucleic_acid<NucleicAcid>::value, "Values are not nucleic acids.");
        using UT = typename std::underlying_type<NucleicAcid>::type;
        return static_cast<NucleicAcid>(~static_cast<UT>(nt) & 0b1111);
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

    template <typename NucleicAcid>
    bool is_certain(NucleicAcid nt){
        return count_ones(nt) == 1;
    }

    DNA char_to_dna[256] = { DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid, 
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
                                 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Gap, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::A,
				 DNA::B, DNA::C, DNA::D,
				 DNA::Invalid, DNA::Invalid, DNA::G,
				 DNA::H, DNA::Invalid, DNA::Invalid,
				 DNA::K, DNA::Invalid, DNA::M,
				 DNA::N, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::R, DNA::S,
                                 DNA::T, DNA::Invalid, DNA::V,
				 DNA::W, DNA::Invalid, DNA::Y,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::A, DNA::B,
				 DNA::C, DNA::D, DNA::Invalid,
				 DNA::Invalid, DNA::G, DNA::H,
				 DNA::Invalid, DNA::Invalid, DNA::K,
				 DNA::Invalid, DNA::M, DNA::N,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::R, DNA::S, DNA::T,
				 DNA::Invalid, DNA::V, DNA::W,
				 DNA::Invalid, DNA::Y, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
                                 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
                                 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
                                 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid, DNA::Invalid, DNA::Invalid,
                                 DNA::Invalid, DNA::Invalid, DNA::Invalid,
				 DNA::Invalid };
    
    char dna_to_char[16] = { '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W',
			     'Y', 'H', 'K', 'D', 'B', 'N' };

    RNA char_to_rna[256] = { RNA::Invalid, RNA::Invalid, RNA::Invalid,
	                     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Gap, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::A,
			     RNA::B, RNA::C, RNA::D,
			     RNA::Invalid, RNA::Invalid, RNA::G,
			     RNA::H, RNA::Invalid, RNA::Invalid,
			     RNA::K, RNA::Invalid, RNA::M,
			     RNA::N, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::R, RNA::S,
			     RNA::Invalid, RNA::U, RNA::V,
			     RNA::W, RNA::Invalid, RNA::Y,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::A, RNA::B,
			     RNA::C, RNA::D, RNA::Invalid,
			     RNA::Invalid, RNA::G, RNA::H,
			     RNA::Invalid, RNA::Invalid, RNA::K,
			     RNA::Invalid, RNA::M, RNA::N,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::R, RNA::S, RNA::Invalid,
			     RNA::U, RNA::V, RNA::W,
			     RNA::Invalid, RNA::Y, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid, RNA::Invalid, RNA::Invalid,
			     RNA::Invalid };

    char rna_to_char[16] = { '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W',
			     'Y', 'H', 'K', 'D', 'B', 'N' };

}
#endif
