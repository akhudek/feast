#ifndef LOCAL_SEGMENTATION_H__
#define LOCAL_SEGMENTATION_H__
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "ereal.h"
#include "zero_order_background.h"
#include "dna_sequence.h"
//#include "similarity_region.h"
#include "nw_model_parameters.h"
#include <btl/logspace.h>

template<class background_model = zero_order_background, typename DT = ereal, typename RDT = ereal >
class local_segmentation {
public:
	typedef DT ftype;

	struct result {
		RDT score;
		int a;
		int b;
		result() : score(0.0), a(-1), b(-1) {}
		result( DT sc, int na, int nb ) : score(sc), a(na), b(nb) {}
	};

protected:
	typedef unsigned int uint;

	// dynamic programming cell, one variable for each state
	struct dpcell{
		DT M;
		DT Ia;
		DT Ib;
	};
	typedef blitz::Array<dpcell,2> dparray;

	// we store several precomputed values at each location
	struct precomputed_values {
		DT T3;
		DT T4;
		DT IpJx2;
		DT IpJm1;
	};

	blitz::Array<precomputed_values,2> pv;
	static int const MAX_EXT = 400;

	// horizontal frontier
	dparray dpr;

	void init_storage();

	// The main DP.
	result maindp( dna_sequence_data &a, dna_sequence_data &b, int i, int j, bool forward, int ext );
	result maindp_filldp( dna_sequence_data &a, dna_sequence_data &b, int i, int j, bool forward, int ext );

	// probability score matrix for DNA sequences in log_sqrt(2) scale
	blitz::Array<DT,2> p;

	// background probabilities
	blitz::TinyVector<DT,4> q;

	// probabilities for gap open and gap extend also in log_sqrt(2) scale
	DT p_a, p_1m2a, p_b, p_1mb, pi_m, pi_i;
	double pr_a, pr_b;

	// where we terminate extensions on masked sequences
	bool stop_on_masked;

	// initialize the scoring model with BLASTN defaults
	void init_pr();

	// compute the stationary distribution
	void compute_sd();

	// initialize a dp matrix location
	void init_location( dparray &dp, int i, int j );

public:
	// default constructor
	local_segmentation();

	// copy constructor
	local_segmentation( local_segmentation const &other );

	// extend forwards
	result extend_forward( dna_sequence &seqa, dna_sequence &seqb, int i, int j, int ext = 200 );

	// extend backwards
	result extend_backward( dna_sequence &seqa, dna_sequence &seqb, int i, int j, int ext = 200 );

	// set all the parameters
	void set_parameters( nw_model_parameters const &mp );

	// set gap open
	void set_popen( double new_popen );

	// set gap extend
	void set_mean_gap_length( double d );

	// set mutation matrix
	void set_p( blitz::Array<double,2> const &newp );

	// set background
	void set_q( blitz::TinyVector<double,4> const &newq );

	// set stop on masked flag
	void set_stop_on_masked( bool );

};


template<class background_model,typename DT, typename RDT>
inline void local_segmentation<background_model,DT,RDT>::set_stop_on_masked( bool m ) {
	stop_on_masked = m;
}

template<class background_model,typename DT, typename RDT>
inline void local_segmentation<background_model,DT,RDT>::set_parameters( nw_model_parameters const &mp ) {
	set_q( mp.q );
	set_p( mp.p );
	set_popen( mp.pr_open );
	set_mean_gap_length( mp.mean_gap_length );
}

template<class background_model, typename DT, typename RDT>
inline void local_segmentation<background_model,DT,RDT>::init_location( dparray &dp, int i, int j ) {
   dp(i,j).M = dp(i,j).Ia = dp(i,j).Ib = 0.0;
}

// initialize storage
template<class background_model, typename DT, typename RDT>
void local_segmentation<background_model,DT,RDT>::init_storage() {
   p.resize(dna_alpha::SIZE,dna_alpha::SIZE);

   // precompute some values
   pv.resize(MAX_EXT+1,MAX_EXT+1);
   dpr.resize(MAX_EXT+1,2);
   for( int i = 1; i <= MAX_EXT; ++i ) {
		for( int j = 1; j <= MAX_EXT; ++j ) {
			pv(i,j).T3 = (double)(i*i + j*j + 2*(i*j - i - j) + 1);
			pv(i,j).T4 = (double)(4*(i*i + j*j) + 8*(i*j - i - j));
			pv(i,j).IpJx2 = (double)(2*(i+j));
			pv(i,j).IpJm1 = (double)(i+j-1);
		}
   }
}

// We assume the BLASTN default scoring scheme with Juke's/Cantor background. All logs are in base sqrt(2).
// Defaults          Adjusted to remove/include the odds portion.
// Match:       1    PrMatch:       log(4^{1}/16)     = -4     remove odds
// Mismatch:   -3    PrMismatch:    log(4^{-3}/16)    = -20    remove odds
// Gap open:   -5    Pr gap open:   log(4^{-5})       = -20    odds cancel, no need to remove
// Gap extend: -2    Pr gap extend: log(4)^{-2})      = -8    odds cancel, no need to remove
// Pr A,T,G,C: 1/4   Pr A,T,G,C:    log(1/4)          = -4     add up odds separately
template<class background_model, typename DT, typename RDT>
void local_segmentation<background_model, DT, RDT>::init_pr() {
   DT ma = 0.25; // match
   DT mm = std::pow( std::sqrt(2.0), -20.0 );
   p =   ma, mm, mm, mm,
	mm, ma, mm, mm,
	mm, mm, ma, mm,
	mm, mm, mm, ma;

   // background probabilities
   q =   ma, ma, ma, ma;

   // precompute transition probabilities
   set_popen( 0.001043033 );
   set_mean_gap_length(1.07);
}

// set mean gap length to n
template<class background_model, typename DT, typename RDT>
void local_segmentation<background_model,DT, RDT>::set_mean_gap_length( double n ) {
   // Since a gap is length 1 by the virtue of opening it, we adjust n down by 1.
   // The remaining gap extension probability corresponds to the mean of a
   // geometric distribution allowing zero extensions.
   assert( n >= 1 );
   pr_b = 1.0-1.0/n;
   p_b = pr_b;
   p_1mb = 1.0 - pr_b;
   compute_sd();
}

// convert to log base sqrt(2)
template<class background_model,typename DT, typename RDT>
void local_segmentation<background_model,DT,RDT>::set_popen( double new_open ) {
   using namespace std;
   assert( new_open <= 1.0 && new_open >= 0 );
   pr_a = new_open;
   p_a = pr_a;
   p_1m2a = 1.0-2.0*pr_a;
   compute_sd();
}

template<class background_model,typename DT, typename RDT>
void local_segmentation<background_model,DT,RDT>::compute_sd() {
   double base = 1.0 - pr_b + 2.0*pr_a;
   pi_m = (1.0 - pr_b)/base;
   pi_i = pr_a/base;

   assert( (double)pi_m >= 0 && (double)pi_m <= 1.0 );
   assert( (double)pi_i >= 0 && (double)pi_i <= 1.0 );
}

// convert each value to log base sqrt(2)
template<class background_model,typename DT, typename RDT>
void local_segmentation<background_model,DT,RDT>::set_q( blitz::TinyVector<double,4> const &newq ) {
   using namespace std;
   for( int i = 0; i < 4; i++ ) {
      q(i) = newq(i);
   }
}

// convert each value to log base sqrt(2)
template<class background_model,typename DT, typename RDT>
void local_segmentation<background_model,DT,RDT>::set_p( blitz::Array<double,2> const &newp ) {
   using namespace std;
   assert( newp.rows() == 4 );
   assert( newp.cols() == 4 );

   for( int i = 0; i < 4; i++ ) for( int j = 0; j < 4; j++ ) p(i,j) = newp(i,j);
}

template<class background_model,typename DT, typename RDT>
typename local_segmentation<background_model,DT,RDT>::result local_segmentation<background_model,DT,RDT>::maindp_filldp( dna_sequence_data &a, dna_sequence_data &b, int si, int sj, bool forward, int ext ) {
	using namespace std;
	int n = (int) a.size(), m = (int)b.size();
	assert( n > 0 );
	assert( m > 0 );
	assert( ext > 0 );
	assert( ext < MAX_EXT );
	int ei = forward ? std::min(si+ext-1,n-1) : std::max(si-ext+1,0);
	int ej = forward ? std::min(sj+ext-1,m-1) : std::max(sj-ext+1,0);

    // indicates which row is the front
    int fr = 1, ba = 0;

    // precompute odds
    closed_interval max_interval_a = forward ? closed_interval(si,ei) : closed_interval(ei,si);
    closed_interval max_interval_b = forward ? closed_interval(sj,ej) : closed_interval(ej,sj);
    background_model oddsA = background_model( q, a, max_interval_a );
    background_model oddsB = background_model( q, b, max_interval_b );

    // full background model
    int last_i = max_interval_a.size() - 1, last_j = max_interval_b.size() - 1;
    DT odds_full = oddsA.pr(closed_interval(0,last_i))*oddsB.pr(closed_interval(0,last_j));

    // declare some constants and variables ahead so we don't keep reallocating them during the loop
	DT odds_a, odds_b;
	int ai = si, bi = sj;

	// initialize top row
	for( int j = 0; j <= ext; ++j ) init_location(dpr,j,ba);

	// we may only start in 1,1 match position
	init_location(dpr,1,fr);
	dpr(1,fr).M =  p((int)a[si],(int)b[sj]);

    // store best score and location in matrix
    DT sc, best_s = 0.0;
    int best_i = 0, best_j = 0;

	// here we skip 1,1 since it is set explicitly
	int ext_i = abs(ei-si)+1, ext_j = abs(ej-sj)+1;
	for( int i = 1, j = 2; i <= ext_i; ++i ) {
		// initialize leftmost column to zero
		init_location(dpr,0,fr);

		for( ; j <= ext_j; ++j ) {
			// Compute odds for remainder of sequence (the random part).
			if( forward ) {
				ai = si+i-1; bi = sj+j-1;
				odds_a = ( i == ei ) ? DT(1.0) : oddsA.pr(closed_interval(i,last_i));
				odds_b = ( j == ej ) ? DT(1.0) : oddsB.pr(closed_interval(j,last_j));
			} else {
				ai = si-i+1; bi = sj-j+1;
				odds_a = ( ai - 1 < ei ) ? DT(1.0) : oddsA.pr(closed_interval(0,last_i - i));
				odds_b = ( bi - 1 < ej ) ? DT(1.0) : oddsB.pr(closed_interval(0,last_j - j));
			}
#ifdef VITERBI_EXTENSIONS
			dpr(j,fr).M =  max(p_1m2a*dpr(j-1,ba).M, p_1mb*max(dpr(j-1,ba).Ia, dpr(j-1,ba).Ib))
					* p((int)a[ai],(int)b[bi]);
			dpr(j,fr).Ia = max(p_a*dpr(j,ba).M, p_b*dpr(j,ba).Ia)*q((int)a[ai]);
			dpr(j,fr).Ib = max(p_a*dpr(j-1,fr).M,p_b*dpr(j-1,fr).Ib)*q((int)b[bi]);
			sc = (dpr(j,fr).Ia + dpr(j,fr).Ib + dpr(j,fr).M)*odds_a*odds_b/odds_full;
#else
			dpr(j,fr).M = (p_1m2a*dpr(j-1,ba).M + p_1mb*(dpr(j-1,ba).Ia + dpr(j-1,ba).Ib))* p((int)a[ai],(int)b[bi]);
			dpr(j,fr).Ia = (p_a*dpr(j,ba).M + p_b*dpr(j,ba).Ia)*q((int)a[ai]);
			dpr(j,fr).Ib = (p_a*dpr(j-1,fr).M + p_b*dpr(j-1,fr).Ib)*q((int)b[bi]);
			sc = (dpr(j,fr).Ia + dpr(j,fr).Ib + dpr(j,fr).M)*odds_a*odds_b/odds_full;
#endif

			// store overall best
			if( sc > best_s ) {
				best_s = sc;
				best_i = i;
				best_j = j;
			}
			//cerr << i << "\t" << j << "\t" << best_s << "\t" << odds_full << endl;

			// stop extension if we encounter repeat
			if( a[ai].masked() && stop_on_masked ) ext_i = i;
			if( b[bi].masked() && stop_on_masked ) ext_j = j;
		} // for j

		j = 1; // reset j

		// swap rows
		fr = 1 - fr; ba = 1 - ba;

   } // for i

   // otherwise backwards
   result my_result(best_s,best_i,best_j);
   return my_result;
}

template<class background_model,typename DT, typename RDT>
typename local_segmentation<background_model,DT,RDT>::result local_segmentation<background_model,DT,RDT>::maindp( dna_sequence_data &a, dna_sequence_data &b, int i, int j, bool forward, int ext ) {
   return maindp_filldp( a, b, i, j, forward, ext );
}

// init pr_a and pr_b so our calcuations work out
template<class background_model,typename DT, typename RDT>
local_segmentation<background_model,DT,RDT>::local_segmentation() : pr_a(0.06), pr_b(0.05), stop_on_masked(false) {
   init_storage();
	init_pr();
}

// init pr_a and pr_b so our calcuations work out
template<class background_model,typename DT, typename RDT>
local_segmentation<background_model,DT,RDT>::local_segmentation( local_segmentation<background_model,DT,RDT> const &o ) : pr_a(o.pr_a), pr_b(o.pr_b), stop_on_masked(false) {
   init_storage();
	p = o.p;
	q = o.q;
	p_a = o.p_a;
	p_1m2a = o.p_1m2a;
	p_b = o.p_b;
	p_1mb = o.p_1mb;
	pi_m = o.pi_m;
	pi_i = o.pi_i;
}


template<class background_model,typename DT, typename RDT>
typename local_segmentation<background_model,DT,RDT>::result local_segmentation<background_model,DT,RDT>::extend_forward( dna_sequence &seqa, dna_sequence &seqb, int i, int j, int ext ) {
   return maindp( seqa.data, seqb.data, i, j, true, ext );
}

template<class background_model,typename DT, typename RDT>
typename local_segmentation<background_model,DT,RDT>::result local_segmentation<background_model,DT,RDT>::extend_backward( dna_sequence &seqa, dna_sequence &seqb, int i, int j, int ext ) {
   return maindp( seqa.data, seqb.data, i, j, false, ext );
}

#endif

