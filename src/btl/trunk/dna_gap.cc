#include "btl/dna_gap.h"

btl::dna_gap_symbol const btl::dna_gap::A('A');
btl::dna_gap_symbol const btl::dna_gap::T('T');
btl::dna_gap_symbol const btl::dna_gap::C('C');
btl::dna_gap_symbol const btl::dna_gap::G('G');
btl::dna_gap_symbol const btl::dna_gap::GAP('-');

boost::array<btl::dna_gap_symbol,5> const btl::dna_gap::alphabet =  { 
btl::dna_gap::A, 
btl::dna_gap::T, 
btl::dna_gap::C, 
btl::dna_gap::G, 
btl::dna_gap::GAP
};

