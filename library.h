#ifndef BIOSYMBOLS_LIBRARY_H
#define BIOSYMBOLS_LIBRARY_H

template<typename T>
struct is_nucleic_acid{
    static const bool value = false;
};

// The full gamut of symbols for the DNA nucleotide alphabet.
enum class DNA : unsigned char { Gap, A, C, M, G, R, S, V, T, W, Y, H, K, D, B, N };

template<>
struct is_nucleic_acid<DNA>{
    static const bool value = true;
};

// The full gamut of symbols for the RNA alphabet.
enum class RNA : unsigned char { Gap, A, C, M, G, R, S, V, U, W, Y, H, K, D, B, N };

template<>
struct is_nucleic_acid<RNA>{
    static const bool value = true;
};

// The full gamut of symbols for the amino acid alphabet.
enum class AA : unsigned char { A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, O, U, B, J, Z, X, Term, Gap };

// Test if `nt` is surely either guanine or cytosine.
bool isGC(DNA nt);

// Test if `nt` is surely a purine.
bool ispurine(DNA nt);

// Test if `nt` is surely a pyrimidine.
bool ispyrimidine(DNA nt);

// Test if `nt` represents an ambiguity symbol.
bool isambiguous(DNA nt);

// Test if `nt` represents a certainly known nucleotide.
bool iscertain(DNA nt);

// Test if `nt` is a symbol representing a gap.
bool isgap(DNA nt);


#endif
