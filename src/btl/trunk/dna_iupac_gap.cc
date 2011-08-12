#include "btl/dna_iupac_gap.h"

btl::dna_iupac_gap::symbol const btl::dna_iupac_gap::alphabet_decode[] 
   = "-GCSTKYBARMVWDHN~";
btl::dna_iupac_gap::symbol const btl::dna_iupac_gap::alphabet[] = {
   btl::dna_iupac_gap::GAP,
   btl::dna_iupac_gap::G,
   btl::dna_iupac_gap::C,
   btl::dna_iupac_gap::S,
   btl::dna_iupac_gap::T,
   btl::dna_iupac_gap::K,
   btl::dna_iupac_gap::Y,
   btl::dna_iupac_gap::B,
   btl::dna_iupac_gap::A,
   btl::dna_iupac_gap::R,
   btl::dna_iupac_gap::M,
   btl::dna_iupac_gap::V,
   btl::dna_iupac_gap::W,
   btl::dna_iupac_gap::D,
   btl::dna_iupac_gap::H,
   btl::dna_iupac_gap::N,
   btl::dna_iupac_gap::end_token };

int const btl::dna_iupac_gap::symbol_ambiguity[] 
   = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };


