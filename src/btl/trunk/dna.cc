#include "btl/dna.h"

btl::dna_symbol const btl::dna::A('A');
btl::dna_symbol const btl::dna::T('T');
btl::dna_symbol const btl::dna::C('C');
btl::dna_symbol const btl::dna::G('G');

boost::array<btl::dna_symbol,4> const btl::dna::alphabet =  { 
   btl::dna::A, 
   btl::dna::T, 
   btl::dna::C, 
   btl::dna::G, 
};

