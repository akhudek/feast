#include "alignment_functions.h"

// append b to a
void append_alignment( pairwise_dna_alignment &align_a, pairwise_dna_alignment &align_b ) {
   append_alignment_sequence( *align_a.a, *align_b.a );
   append_alignment_sequence( *align_a.b, *align_b.b );
   align_a.score += align_b.score;
}

void append_alignment_sequence( dna_alignment_sequence &a, dna_alignment_sequence &b ) {
   for( dna_alignment_sequence_data_citer i = b.data.begin(); i != b.data.end(); i++ ) {
      a.data.push_back( *i );
   }
}

// append gap of length n to a
void append_gap( dna_alignment_sequence &a, int n ) {
   for( int i = 0; i < n; i++ ) a.data.push_back( dna_alignment_alpha::GAP );
}

// append region b to alignment sequence a
void append_region( dna_alignment_sequence &a, dna_sequence_region &b ) {
   for( dna_sequence_region_data_citer i = b.data.begin(); i != b.data.end(); i++ ) {
      a.data.push_back( *i );
   }
}




