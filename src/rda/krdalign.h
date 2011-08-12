#ifndef KRDALIGN_H__
#define KRDALIGN_H__
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <btl/logspace.h>
#include "nw_model_parameters.h"
#include "dna_sequence_region.h"
#include "similarity_region.h"
#include "ereal.h"

class krdalign {
protected:
   typedef unsigned int uint;
   
	// maximum length of L or R block
	unsigned int k, shape_y;
	int shape_z;
	
	// Arrays we use to compute lalign and ralign blocks. We reduce allocation time by reusing these arrays instead of
	// continually allocating new ones.
	typedef blitz::Array< ereal ,2> smatrix;
	smatrix M, I, Sd, Sr; // ,St;
	blitz::Array< ereal ,1> odds_a, odds_b;
	
	void init_storage();
	void init_shape();
	
	int shape_f( int x );
	blitz::Array<int,1> shape_fx, shape_fy;
	
	
	// data structure for main dp, can be packed tighter
	struct dpcell {
		ereal  s;
		uint8_t from;
		int u;
		int l;
	};
	struct dpdata {
		dpcell R;
		dpcell D;
	};
	enum states { no_state = 0,  M_state, Ga_state, Gb_state, R_block, D_block, num_states };
	typedef blitz::Array<dpdata,2> dparray;
	
	void init_cell( dpcell &s );
	void init_location( dparray &dp, int i, int j );
	inline ereal state_score( dparray &dp, int i, int j, uint8_t type );
	
	// The main DP.
	std::pair<similarity_region_deque_ptr,double> maindp( dna_sequence_region_data &a, dna_sequence_region_data &b );
	void maindp_init_borders( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	void maindp_filldp( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	std::pair<similarity_region_deque_ptr,double> maindp_traceback( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b );
	
	// Compute a down diagonal block. Just swap a and b to compute right diagonal blocks. 
	void hmmdp( dna_sequence_region_data &a, uint ai, uint aj, dna_sequence_region_data &b, uint bi, uint bj, smatrix &S, int at );
	
	enum at_values { BETWEEN_BLOCKS = 0, START = 1, END = 2 };
	
	// probability score matrix for DNA sequences in log_sqrt(2) scale
	smatrix p;
	
	// background probabilities
	blitz::TinyVector< ereal ,4> q;
	
	// probabilities for state transitions
	ereal p_a, p_b, p_1m2a, p_1mb;
	
	// store original parameters for accurate computation of the stationary distribution
	double original_a, original_b;
	
	// initialize the scoring model with BLASTN defaults
	void init_pr();
	
	// our starting probabilities
	double pr_red, pr_m, pr_i;
	
	// precompute starting probabilities
	void compute_sd();
	
	
public:
	// default constructor
	krdalign();
	
	std::pair<similarity_region_deque_ptr,double> align( dna_sequence_region &seqa, dna_sequence_region &seqb );
	
	// set k
	void set_k( unsigned int k );
	
	// Set all parameters.
	void set_parameters( nw_model_parameters &mp );
	
	// set gap open (Pr[open gap of length >=1 in one sequence]
	void set_popen( double new_popen );
	
	// set gap extend
	void set_mean_gap_length( double n );
	
	// set mutation matrix
	void set_p( blitz::Array<double,2> &newp );
	
	// set background 
	void set_q( blitz::TinyVector<double,4> &newq );
};

// Inline functions follow.
// --------------------------------------------

inline void krdalign::init_cell( dpcell &c ) {
   c.s = 0.0;
   c.from = no_state;
   c.u = -1;
   c.l = -1;
}

inline void krdalign::init_location( dparray &dp, int i, int j ) {
   init_cell( dp(i,j).R );
   init_cell( dp(i,j).D );
}

inline ereal krdalign::state_score( dparray &dp, int i, int j, uint8_t type ) {
   assert( type < num_states );
   assert( type != no_state );
   switch( type ) {
      case R_block:
         return dp(i,j).R.s;
   }
   return dp(i,j).D.s;
}

inline int krdalign::shape_f( int x ) {
   using namespace std;
   double dk = (double)k;
   double dz = (double)shape_z;
   double dy = (double)shape_y;
   double dx = (double) x;
   return (int) (2.0 + (dk-2.0)*pow(dz-dx, dy)/pow(dz,dy));
}

#endif

