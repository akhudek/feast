#ifndef __PAIRWISE_DNA_ALIGNMENT_H_
#define __PAIRWISE_DNA_ALIGNMENT_H_
#include "dna_alignment_sequence.h"

struct pairwise_dna_alignment {
   dna_alignment_sequence_ptr a;
   dna_alignment_sequence_ptr b;
   double score;

   pairwise_dna_alignment( dna_alignment_sequence_ptr na, dna_alignment_sequence_ptr nb, double s ) :
      a(na),
      b(nb),
      score(s) {}

   pairwise_dna_alignment() : score(0.0) {}
};

#endif
