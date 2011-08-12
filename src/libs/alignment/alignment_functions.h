#ifndef ALIGNMENT_FUNCTIONS_H_
#define ALIGNMENT_FUNCTIONS_H_

#include "dna_alignment_sequence.h"
#include "dna_sequence_region.h"
#include "pairwise_dna_alignment.h"

void append_alignment( pairwise_dna_alignment &a, pairwise_dna_alignment &b );

void append_gap( dna_alignment_sequence &a, int i );

void append_region( dna_alignment_sequence &a, dna_sequence_region &b );

void append_alignment_sequence( dna_alignment_sequence &a, dna_alignment_sequence &b );

#endif
