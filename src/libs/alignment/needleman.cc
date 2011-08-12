#include "needleman.h"

#include "ereal.h"
#include <cassert>
#include <algorithm>
#include <blitz/tinyvec-et.h>
#include <btl/logspace.h>

// We assume the BLASTN default scoring scheme with Juke's/Cantor background. All logs are in base sqrt(2).
// Defaults          Adjusted to remove/include the odds portion.
// Match:       1   4
// Mismatch:   -3   -12
// Gap open:   -5   -20
// Gap extend: -2   -8
void needleman::init_pr() {
   // match/mismatch probabilities
   // pno_open should be rolled in here, but it's value is -0.006 which is
   // barely representable by the precisian we use here. Ignore for now.
   s =   4.0, 0.015625, 0.015625, 0.015625,
         0.015625, 4.0, 0.015625, 0.015625,
         0.015625, 0.015625, 4.0, 0.015625,
         0.015625, 0.015625, 0.015625, 4.0;

   // gap open
   popen = 0.001043033;

   // gap extend
   pext = 0.0625;
   pno_ext = 1.0-0.0625;
}

void needleman::set_parameters( nw_model_parameters const &mp ){
	set_mean_gap_length( mp.mean_gap_length );
	set_s( mp.pr_open, mp.p, mp.q );
}

// set mean gap length to n
void needleman::set_mean_gap_length( double n ) {
   // Since a gap is length 1 by the virtue of opening it, we adjust n down by 1.
   // The remaining gap extension probability corresponds to the mean of a
   // geometric distribution allowing zero extensions.
   pext =  1.0-1.0/n;
   pno_ext =  1.0/n;
	}

// computes score matrix from match probabilities and background
void needleman::set_s( double new_popen, blitz::Array<double,2> const &p, blitz::TinyVector<double,4> const &q ) {
   using namespace std;
   assert( p.rows() == 4 );
   assert( p.cols() == 4 );

   popen = new_popen;
	double pno_open = 1.0-2*new_popen;
   for( int i = 0; i < 4; i++ ) for( int j = 0; j < 4; j++ ) s(i,j) = pno_open*p(i,j)/q(i)/q(j);
}

void needleman::maindp_init_borders( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   int n = (int) a.size(), m = (int)b.size();
   assert( (n>=1) && (m>=1) );

   // top corner is the only valid starting point
   init_cell(dp(0,0).M); init_cell(dp(0,0).Ga); init_cell(dp(0,0).Gb);
   dp(0,0).M.s = 1.0;

   // Initialize borders. Scores represent opening a big gap.
   init_cell(dp(1,0).M); init_cell(dp(1,0).Ga); init_cell(dp(1,0).Gb);
   dp(1,0).Gb.from = M_state; dp(1,0).Gb.s = popen*pno_ext;
   for( int i = 2; i <= n; i++ ) {
      init_cell(dp(i,0).M); init_cell(dp(i,0).Ga); init_cell(dp(i,0).Gb);
      dp(i,0).Gb.from = Gb_state; dp(i,0).Gb.s = dp(i-1,0).Gb.s*pext;
   }
   init_cell(dp(0,1).M); init_cell(dp(0,1).Ga); init_cell(dp(0,1).Gb);
   dp(0,1).Ga.from = M_state; dp(0,1).Ga.s = popen*pno_ext;
   for( int j = 2; j <= m; j++ ) {
      init_cell(dp(0,j).M); init_cell(dp(0,j).Ga); init_cell(dp(0,j).Gb);
      dp(0,j).Ga.from = Ga_state; dp(0,j).Ga.s = dp(0,j-1).Ga.s*pext;
   }
}

void needleman::maindp_filldp( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   int n = (int) a.size(), m = (int)b.size();
	ereal score, match_cost;
   for( int i = 1; i <= n; i++ ) {
   for( int j = 1; j <= m; j++ ) {
      // first consider a match, default is come from match
      match_cost = s((int)a[(uint)i-1],(int)b[(uint)j-1]);
      dp(i,j).M.s = dp(i-1,j-1).M.s*match_cost;
      dp(i,j).M.from = M_state;

      // consider coming from a gap in A
      score = dp(i-1,j-1).Ga.s*match_cost;
      if( score > dp(i,j).M.s ) {
         dp(i,j).M.s = score;
         dp(i,j).M.from = Ga_state;
      }

      // consider coming from a gap in B
      score = dp(i-1,j-1).Gb.s*match_cost;
      if( score > dp(i,j).M.s ) {
         dp(i,j).M.s = score;
         dp(i,j).M.from = Gb_state;
      }

      // now consider extending or opening a gap in A, default is open a gap
      dp(i,j).Ga.s = dp(i,j-1).M.s*popen*pno_ext;
      dp(i,j).Ga.from = M_state;

      // extend a gap
      score = dp(i,j-1).Ga.s*pext;
      if( score > dp(i,j).Ga.s ) {
         dp(i,j).Ga.s = score;
         dp(i,j).Ga.from = Ga_state;
      }

      // now consider extending or opening a gap in B, default is open a gap
      dp(i,j).Gb.s = dp(i-1,j).M.s*popen*pno_ext;
      dp(i,j).Gb.from = M_state;

      // extend a gap
      score = dp(i-1,j).Gb.s*pext;
      if( score > dp(i,j).Gb.s ) {
         dp(i,j).Gb.s = score;
         dp(i,j).Gb.from = Gb_state;
      }
   } // for j
   } // for i
}

pairwise_dna_alignment needleman::maindp_traceback( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   using namespace std;
   int n = (int) a.size(), m = (int)b.size();

   // alignment sequences
   dna_alignment_sequence_ptr newa = new_dna_alignment_sequence(), newb = new_dna_alignment_sequence();

   // determine starting block
   int i = n, j = m, next_i = n, next_j = m;
   uint8_t type = M_state, next_type = no_state;
   if( state_score(dp,i,j,type) < dp(i,j).Ga.s ) type = Ga_state;
   if( state_score(dp,i,j,type) < dp(i,j).Gb.s ) type = Gb_state;
   ereal alignment_score = state_score(dp,i,j,type);

   while( (i > 0) || (j > 0) ) {
      assert( type != no_state );
      assert( type < num_states );
      switch( type ) {
         case M_state:
            next_i = i - 1;
            next_j = j - 1;
            next_type = dp(i,j).M.from;
            newa->data.push_front( a[(uint)i-1] );
            newb->data.push_front( b[(uint)j-1] );
            break;
         case Ga_state:
            next_i = i;
            next_j = j - 1;
            next_type = dp(i,j).Ga.from;
            newa->data.push_front( dna_alignment_alpha::GAP );
            newb->data.push_front( b[(uint)j-1] );
            break;
         case Gb_state:
            next_i = i - 1;
            next_j = j;
            next_type = dp(i,j).Gb.from;
            newa->data.push_front( a[(uint)i-1] );
            newb->data.push_front( dna_alignment_alpha::GAP );
      }
      i = next_i;
      j = next_j;
      type = next_type;
   }
   return pairwise_dna_alignment( newa, newb, alignment_score.as_base() );
}

pairwise_dna_alignment needleman::maindp( dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   int n = (int) a.size(), m = (int)b.size();
   dparray dp(n+1,m+1);
   maindp_init_borders( dp, a, b );
   maindp_filldp( dp, a, b );
   return maindp_traceback( dp, a, b );
}

// only need to initialize score scheme
needleman::needleman() {
   s.resize(dna_alpha::SIZE,dna_alpha::SIZE);
   init_pr();
}

pairwise_dna_alignment needleman::align( dna_sequence_region &seqa, dna_sequence_region &seqb ) {
   return maindp( seqa.data, seqb.data );
}

std::pair<nw_model_parameters_ptr,ereal > needleman::estimate( dna_sequence_region &seqa, dna_sequence_region &seqb ) {
	using namespace std;
	pairwise_dna_alignment aln = align(seqa,seqb);
	counts cnt;
	do_statistics(aln,cnt);
	nw_model_parameters_ptr nwp = update_probabilities(cnt);
	return make_pair(nwp,aln.score);
}

void needleman::do_statistics( pairwise_dna_alignment &align, needleman::counts &cnt ) {
	assert( align.a->data.size() == align.b->data.size() );

	int gap_a = -1, gap_b = -1;
	for( int c = 0; c < (int)align.a->data.size(); ++c ) {
		dna_alignment_alpha::symbol c_a = align.a->data[c], c_b = align.b->data[c];

		// found a pair!
		if( c_a != dna_alignment_alpha::GAP && c_b != dna_alignment_alpha::GAP ) {
			int base_a = dna_alpha::symbol(c_a).index();
			int base_b = dna_alpha::symbol(c_b).index();
			cnt.pairs( base_a, base_b) += 1;
			cnt.bases_a(base_a) += 1;
			cnt.bases_b(base_b) += 1;

			// closing a gap in A
			if( gap_a >= 0 ) {
				cnt.gap_lengths_a += c - gap_a + 1;
				gap_a = -1;
			}

			// closing a gap in B
			if( gap_b >= 0 ) {
				cnt.gap_lengths_b += c - gap_b + 1;
				gap_b = -1;
			}

			// found a gap in A
		} else if( c_a == dna_alignment_alpha::GAP && c_b != dna_alignment_alpha::GAP ) {
			cnt.bases_b((int)dna_alpha::symbol(c_b)) += 1;

			// gap is not open
			if( gap_a == -1 ) {
				gap_a = c;
				cnt.gaps_a += 1;
			} // otherwise nothing to do

			// closing a gap in B
			if( gap_b >= 0 ) {
				cnt.gap_lengths_b += c - gap_b + 1;
				gap_b = -1;
			}
			// found a gap in B
		} else if( c_a != dna_alignment_alpha::GAP && c_b == dna_alignment_alpha::GAP ) {
			cnt.bases_a((int)dna_alpha::symbol(c_a)) += 1;

			// closing a gap in A
			if( gap_a >= 0 ) {
				cnt.gap_lengths_a += c - gap_a + 1;
				gap_a = -1;
			}

			// gap is not open
			if( gap_b == -1 ) {
				gap_b = c;
				cnt.gaps_b += 1;
			} // otherwise nothing to do
		} else {
			using namespace std;
			cerr << c << "\t" << c_a << "\t" << c_b << endl;
			assert( 0 && "corrupt alignment" );
		}

	} // for c

	// finished region, count any last gaps
	// closing a gap in A
	if( gap_a >= 0 ) {
		cnt.gap_lengths_a += align.a->data.size() - gap_a;
		gap_a = -1;
	}

	// closing a gap in B
	if( gap_b >= 0 ) {
		cnt.gap_lengths_b += align.a->data.size() - gap_b;
		gap_b = -1;
	}

}

nw_model_parameters_ptr needleman::update_probabilities( needleman::counts &cnt ) {
	using namespace std;
	using namespace blitz;

	// get combined base count for both sequences
	base_count all_bases;
	all_bases = cnt.bases_a + cnt.bases_b;
	int total_bases = sum( all_bases );

	nw_model_parameters_ptr newp( new nw_model_parameters() );

	// compute p
	int total_pairs = sum( cnt.pairs(Range::all(),Range::all()) );
	for( int i = 0; i < dna_alpha::SIZE; ++i ) {
		for( int j = 0; j < dna_alpha::SIZE; ++j ) {
			newp->p(i,j) = (double) cnt.pairs(i,j)/(double)total_pairs;
		} // j
	} // i

	for( int i = 0; i < dna_alpha::SIZE; ++i ) {
		newp->q(i) = (double) all_bases(i)/(double)total_bases;
	}

	// gap open
	int total_gaps = cnt.gaps_a + cnt.gaps_b;
	newp->pr_open = (double)total_gaps/(double)total_bases;

	// mean gap length
	newp->mean_gap_length = (double)(cnt.gap_lengths_a+cnt.gap_lengths_b)/(double)total_gaps;
	return newp;
}

