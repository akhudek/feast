#include "sequence_statistics_functions.h"
#include "dna_sequence.h"

double fraction_repeats( sequence_statistics const &stat, closed_interval const &f ) {
   assert( f.a >= 0 && f.a < (int)stat.repeat_bases.size() );
   assert( f.b >= 0 && f.b < (int)stat.repeat_bases.size() );
   assert( f.b >= f.a );
   unsigned int repeats = ( f.a > 0 ) ? stat.repeat_bases[f.b] - stat.repeat_bases[f.a-1] : stat.repeat_bases[f.b];
   return (double)repeats/(double)f.size();
}

blitz::TinyVector<unsigned int,4> composition_frequency( sequence_statistics const &s, closed_interval const &f ) {
   assert( f.a >= 0 && f.a < (int)s.repeat_bases.size() );
   assert( f.b >= 0 && f.b < (int)s.repeat_bases.size() );
	assert( f.b >= f.a );
   blitz::TinyVector<unsigned int,4> q;
   if( f.a > 0 ) {
      q(dna_alpha::A) =  s.A[f.b] - s.A[f.a-1];
      q(dna_alpha::T) =  s.T[f.b] - s.T[f.a-1];
      q(dna_alpha::C) =  s.C[f.b] - s.C[f.a-1];
      q(dna_alpha::G) =  s.G[f.b] - s.G[f.a-1];
   } else {
      q(dna_alpha::A) = s.A[f.b];
      q(dna_alpha::T) = s.T[f.b];
      q(dna_alpha::C) = s.C[f.b];
      q(dna_alpha::G) = s.G[f.b];
   }
   return q;
}

