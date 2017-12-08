//
// Created by Ben Ward on 25/09/2017.
//

#include "biosymbols.h"
#include <iostream>

using namespace biosymbols;

int main(){
    std::cout << (DNA::A == DNA::A) << std::endl;
    std::cout << is_certain(DNA::A) << std::endl;

    using UT = typename std::underlying_type<DNA>::type;

    char letter = 'T';
    DNA nuc = char_to_dna[letter];
    char outletter = dna_to_char[(UT)nuc];

    std::cout << "The letter " << letter << " is the nucleotide " << outletter << std::endl;
    return 0;
}
