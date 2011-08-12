#ifndef SEQUENCE_STATISTICS_FUNCTIONS_H__
#define SEQUENCE_STATISTICS_FUNCTIONS_H__
#include "sequence_statistics.h"
#include "closed_interval.h"
#include <blitz/tinyvec.h>

template< typename Sequence>
sequence_statistics_ptr measure_statistics( Sequence const &seq ) {
   typedef typename Sequence::data_type seq_data;
   typedef typename Sequence::alphabet alpha;

   sequence_statistics_ptr stats( new sequence_statistics() );
   unsigned int count = 0, a = 0, t = 0, c = 0, g = 0;
   for( typename seq_data::const_iterator i = seq.data.begin(); i != seq.data.end(); ++i ) {
      if( i->masked() ) ++count;
      if( *i == alpha::A ) ++a;
      if( *i == alpha::T ) ++t;
      if( *i == alpha::C ) ++c;
      if( *i == alpha::G ) ++g;

      stats->repeat_bases.push_back( count );
      stats->A.push_back( a );
      stats->T.push_back( t );
      stats->C.push_back( c );
      stats->G.push_back( g );
   }

   // measure triples
   /*
	assert( seq.data.size() > 2 );
   typename alpha::symbol last2 = seq.data[0], last1 = seq.data[1];
   triplet_table<unsigned int> tri_count;
   stats->tri.push_back(tri_count);
   stats->tri.push_back(tri_count);
   for( typename seq_data::const_iterator i = seq.data.begin()+2; i != seq.data.end(); ++i ) {
      tri_count(*i, last1, last2) += 1;
      stats->tri.push_back(tri_count);

      last2 = last1;
      last1 = *i;
   }
	 */ // leave out for now since it takes up a huge amount of memory!
   return stats;
}

double fraction_repeats( sequence_statistics const &stats, closed_interval const &f );
blitz::TinyVector<unsigned int,4> composition_frequency( sequence_statistics const &s, closed_interval const &f );


#endif
