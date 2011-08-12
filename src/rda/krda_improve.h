#ifndef KRDA_IMPROVE_H__
#define KRDA_IMPROVE_H__
// Performs RDA followed by specified alignment within block.

#include "krdalign.h"

#include <boost/tuple/tuple.hpp>
#include "alignment_functions.h"

template< typename Aligner >
class krda_improve {
protected:
	krdalign rda;
	Aligner sa;
	bool block_markers;
	
public:
	// Default constructor.
	krda_improve();
	
	// Set k.
	void set_k( unsigned int k );
	
	// Flag indicating whether we preserve block markers.
	void set_block_marker( bool bm );
	
	// Align two sequences.
	pairwise_dna_alignment align( dna_sequence_region &seqa, dna_sequence_region &seqb );
	
	// Set all parameters.
	void set_parameters( nw_model_parameters &mp );
	
	// Set mean gap length.
	void set_mean_gap_length( double n );
	
	// Set popen, p, and q.
	void set_s( double popen, blitz::Array<double,2> &p, blitz::TinyVector<double,4> &q );
};

// Constructor.
template< typename Aligner >
krda_improve<Aligner>::krda_improve() : block_markers(false) {}

template< typename Aligner >
void krda_improve<Aligner>::set_parameters( nw_model_parameters &mp ) {
	rda.set_parameters( mp );
	sa.set_parameters( mp );
}

template< typename Aligner >
void krda_improve<Aligner>::set_mean_gap_length( double n ) {
   rda.set_mean_gap_length( n );
   sa.set_mean_gap_length( n );
}

template< typename Aligner >
void krda_improve<Aligner>::set_s( double popen, blitz::Array<double,2> &p, blitz::TinyVector<double,4> &q ) {
   sa.set_s( popen, p, q );
   rda.set_popen( popen );
   rda.set_p( p );
   rda.set_q( q );
}

template< typename Aligner >
void krda_improve<Aligner>::set_k( unsigned int k ) {
   rda.set_k( k );
}

template< typename Aligner >
void krda_improve<Aligner>::set_block_marker( bool bm ) {
   block_markers = bm;
}

template< typename Aligner >
pairwise_dna_alignment krda_improve<Aligner>::align( dna_sequence_region &seqa, dna_sequence_region &seqb ) {
   // First perform RDA
   double rda_score = 0.0;
	
   similarity_region_deque_ptr rdregions;
   boost::tie(rdregions,rda_score) = rda.align( seqa, seqb );
   
   // Build final alignment by aligning each identified region with Needleman-Wunsch.
   dna_alignment_sequence_ptr alignmenta = new_dna_alignment_sequence();
   dna_alignment_sequence_ptr alignmentb = new_dna_alignment_sequence();
   pairwise_dna_alignment final = pairwise_dna_alignment( alignmenta, alignmentb, 0.0 );
	
   for( similarity_region_deque_iter i = rdregions->begin(); i != rdregions->end(); i++ ) {
      // gap in first sequence
      if( i->first.a == -1 ) {
         assert( i->second.a != -1 );
         append_gap( *final.a, i->second.size() );
         dna_sequence_region bsubseq( seqb, i->second );
         append_region( *final.b, bsubseq );
			
		// gap in second sequence
      } else if( i->second.a == -1 ) {
         assert( i->first.a != -1 );
         append_gap( *final.b, i->first.size() );
         dna_sequence_region asubseq( seqa, i->first );
         append_region( *final.a, asubseq );
			
		// do alignment
      } else {
         dna_sequence_region asubseq(seqa, i->first ), bsubseq(seqb, i->second );
         pairwise_dna_alignment result = sa.align( asubseq, bsubseq );
         append_alignment( final, result );
      }
		
      if( block_markers ) { 
         // add spacer so we know where segments end
         final.a->data.push_back( dna_alignment_alpha::GAP );
         final.b->data.push_back( dna_alignment_alpha::GAP );
      }
   }
   final.score = rda_score;
   return final;
}


#endif
