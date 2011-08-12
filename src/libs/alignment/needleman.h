#ifndef NEEDLEMAN_H__
#define NEEDLEMAN_H__
#include "ereal.h"
#include <blitz/array.h>
#include <btl/logspace.h>
#include "dna_sequence_region.h"
#include "pairwise_dna_alignment.h"
#include "nw_model_parameters.h"

class needleman {
protected:
	typedef unsigned int uint;

	// data structure for main dp
	struct dpcell {
		ereal s;
		uint8_t from;
	};
	struct dpdata {
		dpcell M;     // score for match state
		dpcell Ga;    // score for gap in A
		dpcell Gb;    // score for gap in B
	};
	enum states { no_state = 0, M_state, Ga_state, Gb_state, num_states };
	typedef blitz::Array<dpdata,2> dparray;

	// The main DP.
	pairwise_dna_alignment maindp( dna_sequence_region_data &a, dna_sequence_region_data &b );
	void  maindp_init_borders( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	void  maindp_filldp( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	pairwise_dna_alignment maindp_traceback( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	void init_cell( dpcell &c );
	ereal state_score( dparray &dp, int i, int j, uint8_t state );

	// score matrix for DNA sequences in log_sqrt(2) scale
	blitz::Array<ereal,2> s;

	// log-odds score for gap open and gap extend also in log_sqrt(2) scale (pno_open taken into score matrix)
	ereal popen, pext, pno_ext;

	// initialize the scoring model with BLASTN defaults
	void init_pr();

	// for training
	typedef blitz::Array<int,2> pair_count;
	typedef blitz::TinyVector<int,dna_alpha::SIZE> base_count;

	struct counts {
		pair_count pairs;
		int gaps_a;
		int gaps_b;
		int gap_lengths_a;
		int gap_lengths_b;

		base_count bases_a;
		base_count bases_b;

		counts() : gaps_a(0), gaps_b(0), gap_lengths_a(0), gap_lengths_b(0) {
			// adjust for other models
			pairs.resize(dna_alpha::SIZE,dna_alpha::SIZE);
			pairs = 0;
			bases_a = 0;
			bases_b = 0;
		}
	};

	void do_statistics( pairwise_dna_alignment &align, counts &cnt );
	nw_model_parameters_ptr update_probabilities( counts &cnt );

public:
	// default constructor
	needleman();

	pairwise_dna_alignment align( dna_sequence_region &seqa, dna_sequence_region &seqb );

	void set_parameters( nw_model_parameters const &mp );

	void set_mean_gap_length( double n );

	void set_s( double popen, blitz::Array<double,2> const &p, blitz::TinyVector<double,4> const &q );

	std::pair<nw_model_parameters_ptr,ereal > estimate( dna_sequence_region &seqa, dna_sequence_region &seqb );
};

inline void needleman::init_cell( dpcell &c ) {
   c.s = 0.0;
   c.from = no_state;
}

inline ereal needleman::state_score( dparray &dp, int i, int j, uint8_t state ) {
   assert( state != no_state );
   assert( state < num_states );
   switch( state ) {
      case M_state:
         return dp(i,j).M.s;
      case Ga_state:
         return dp(i,j).Ga.s;
   }
   return dp(i,j).Gb.s;
}

#endif

