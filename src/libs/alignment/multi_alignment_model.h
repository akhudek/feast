#ifndef MULTI_ALIGNMENT_MODEL_H__
#define MULTI_ALIGNMENT_MODEL_H__

#include "pair_hmm.h"
#include "dna_alignment_sequence.h"
#include "dna_sequence_region.h"
#include "nw_model_parameters.h"
#include "zero_order_background.h"
#include "fragment.h"
#include "extension_store.h"
#include "reverse_dna_sequence_region.h"
#include "algorithm_extra.h"
#include <btl/logspace.h>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <list>

#include "point.h"
#include "diagonal_store.h"

class multi_alignment_model {
protected:
	typedef unsigned int uint;

	int nmodels;

	typedef pair_hmm< dna_alpha, ereal > APHMM;
	APHMM PHMM;

	enum states { state_switch = 0, kstates };
	enum sub_stats { state_S = 0, state_M, state_Ia, state_Ib, state_E, SUBSIZE };

	extension_store store;
	diagonal_store dstore;

	blitz::TinyVector<double,4> q;

	bool dirty;

	bool stop_at_repeat;

public:
	multi_alignment_model( int i );

	void start_training();

	ereal end_training();

	template<typename ITERATOR>
	ereal constrained_bw_train( ITERATOR i, ITERATOR end, int nthreads );

	template<typename D1, typename D2, typename O1, typename O2, typename OUTPUT_ITER>
	ereal extend_from_point( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, point p, OUTPUT_ITER oit, double drop_threshold, double cut_threshold );

	template<typename D1, typename D2, typename O1, typename O2>
	ereal ungapped_extend_from_point( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, point p, double drop_threshold );

	int num_models() const;

	void disable_random_model();

	nw_model_parameters_ptr get_parameters( int k ) const;

	void set_parameters( int k, nw_model_parameters const &p, int t );

	void set_default_parameters( blitz::TinyVector<double,4> q );

	void set_stop_at_repeat( bool t );

	template< typename DIST >
	void set_region_dist( DIST const &rdist );

	double get_region_pr( int k ) const;

	double get_region_length( int k ) const;

	void set_random_distribution( blitz::TinyVector<double,4> q );
	blitz::TinyVector<double,4> get_random_distribution() const;

	void set_initial_background( blitz::TinyVector<double,4> q );

	void set_cut_threshold( double t );

	void pretty_print_parameters( std::ostream &out ) const;

	void write_parameters( std::ostream &out ) const;

	void read_parameters( std::istream &in );

	void reset();
};


template<typename ITER>
ereal multi_alignment_model::constrained_bw_train( ITER start, ITER end, int nthreads = 1 ) {
	if( dirty ) { PHMM.optimize(); dirty = false; }
	return PHMM.constrained_train_bw_mt<dna_sequence_region,dna_sequence_region>(nthreads, start,end,40,true);
}


/** Extend in two directions from a point.
 *  - D1 and D2 are expected to be sequence pointers.
 *  - O1 and O2 are odds models for sequences D1 and D2 respectively.
 *  - Point p is the point to start the extension from. The x coord corresponds to D1 and
 *    y to D2.
 *  - Anchor points are sent to oit (for now, maybe using the extensions themselves is better
 *    later on).
 */
template<typename D1, typename D2, typename O1, typename O2, typename OUTPUT_ITER>
ereal multi_alignment_model::extend_from_point( D1 &seqA, O1 &oddsA, D2 &seqB, O2 &oddsB, point p, OUTPUT_ITER oit, double drop_threshold, double cut_threshold ) {
	using namespace std;
	using namespace boost;
	if( dirty ) { PHMM.optimize(); dirty = false; }
	try {
		int max_pos = -1;
		ereal rsc = 1.0;
		extension_store::extension_ptr ext( new extension_store::extension() );

		// Do reverse extension.
		if( p.x > 0 && p.y > 0 ) {
			point rp( p.x-1, p.y-1 );
			extension_store::bounds rbounds = store.get_bounds( rp, false );
			reverse_dna_sequence_region_ptr regAr( new reverse_dna_sequence_region(seqA, closed_interval(0,rp.x)));
			reverse_dna_sequence_region_ptr regBr( new reverse_dna_sequence_region(seqB, closed_interval(0,rp.y)));
			if( stop_at_repeat ) { truncate_at_repeat(*regAr); truncate_at_repeat(*regBr); }

			vector<extension_store::extension::row> reverse_rows;
			tie( max_pos, rsc )= PHMM.forward_extension( *regAr, oddsA, *regBr, oddsB, drop_threshold, back_inserter(reverse_rows), rbounds );

			// Copy data to forward extension and set top range. We also chop off the extension tail.
			if( max_pos >= 0 ) {
			for( vector<extension_store::extension::row>::reverse_iterator i = reverse_rows.rbegin() + (reverse_rows.size() - max_pos - 1);
			     i != reverse_rows.rend(); ++i ) ext->fill.push_back(*i);
			}

		}

		unsigned int reverse_points = ext->fill.size();

		if( ext->fill.size() > 0 ) ext->range.a = p.y - reverse_points;
		else ext->range.a = p.y;

		// Do forward extension.
		extension_store::bounds fbounds = store.get_bounds( p, true );
		dna_sequence_region_ptr regA(new dna_sequence_region(seqA, closed_interval(p.x,seqA->data.size()-1) ));
		dna_sequence_region_ptr regB(new dna_sequence_region(seqB, closed_interval(p.y,seqB->data.size()-1) ));
		if( stop_at_repeat ) { truncate_at_repeat(*regA); truncate_at_repeat(*regB); }


		ereal fsc = 1.0;

		tie(max_pos, fsc) = PHMM.forward_extension( *regA, oddsA, *regB, oddsB, drop_threshold, back_inserter(ext->fill), fbounds );

		// If we have a forward extension, we chop the tail.
		if( max_pos >= 0 ) ext->fill.resize(reverse_points + max_pos+1);

		if( ext->fill.size() > 0 ) ext->range.b = p.y + ext->fill.size() - reverse_points - 1;
		else return 0.0;

		// Ignore extensions that score poorly.
		if( (fsc*rsc).as_base() < cut_threshold ) {
#ifdef DEBUG_OUTPUT
			cerr << "extension too poor" << endl;
#endif
			return 0.0;
		}

		// If we have no left extension, then anchor point is left end.
		int startA = ( reverse_points == 0 ) ? ext->fill.front().range.a : ext->fill.front().best_col ;

		// If we have no right extension, then anchor point is right end.
		int endA = (reverse_points == ext->fill.size() ) ? ext->fill.back().range.b : ext->fill.back().best_col;

		// Get anchor points from extension.
		if( ext->fill.size() <= 200 ) {
			fragment_ptr f( new  fragment() );
			f->add( closed_interval( startA, endA ) );
			f->add( ext->range );
			*oit = f;
			++oit;
		} else {
			closed_interval intA(-1, startA-1), intB(-1 ,ext->range.a-1);

			for( unsigned int i = 1; i < ext->fill.size()-99; ++i ) {
				// We put in an extra check to make sure we've made progress on x.
				if( intA.b+100 < ext->fill[i].best_col &&  (endA - ext->fill[i].best_col) > 100  ) {
					intA.a = intA.b+1; intA.b = ext->fill[i].best_col;
					intB.a = intB.b+1; intB.b = ext->range.a+i;
					fragment_ptr f( new  fragment() );
					f->add( intA );
					f->add( intB );
					*oit = f;
					++oit;
				}
			}
			// Add last segment if it makes progress.
			if( intA.b+1 < endA ) {
				intA.a = intA.b+1; intA.b = endA;
				intB.a = intB.b+1; intB.b = ext->range.b;
				fragment_ptr f( new  fragment() );
				f->add( intA );
				f->add( intB );
				*oit = f;
				++oit;
			}
		}

		// Add extension to store.
		store.add( ext );

		return (fsc*rsc);

	} catch( point_covered_exception &e ) {
#ifdef DEBUG_OUTPUT
		cerr << "detected_overlap" << endl;
#endif
		return 0.0;
	}
}


/** Extend in two directions from a point.
 *  - D1 and D2 are expected to be sequence pointers.
 *  - O1 and O2 are odds models for sequences D1 and D2 respectively.
 *  - Point p is the point to start the extension from. The x coord corresponds to D1 and
 *    y to D2.
 *  - Anchor points are sent to oit (for now, maybe using the extensions themselves is better
 *    later on).
 */
template<typename D1, typename D2, typename O1, typename O2>
ereal multi_alignment_model::ungapped_extend_from_point( D1 &seqA, O1 &oddsA, D2 &seqB, O2 &oddsB, point p , double drop_threshold ) {
	using namespace std;
	using namespace boost;
	using namespace btl;
	if( dirty ) { PHMM.optimize(); dirty = false; }
	if( dstore.overlaps_diagonal(p) ) {
		//cerr << p << " overlaps" << endl;
		return 0.0;
	}
	int dist = 0, distf = 0;
	ereal max_val = 1.0, max_valf = 1.0;

	// Do reverse extension.
	if( p.x > 0 && p.y > 0 ) {
		point rp( p.x-1, p.y-1 );
		closed_interval rrow_interval(0,rp.x);
		reverse_dna_sequence_region_ptr regAr( new reverse_dna_sequence_region(seqA, rrow_interval));
		reverse_dna_sequence_region_ptr regBr( new reverse_dna_sequence_region(seqB, closed_interval(0,rp.y)));
		if( stop_at_repeat ) { truncate_at_repeat(*regAr); truncate_at_repeat(*regBr); }
		tie( max_val, dist ) = PHMM.ungapped_forward_extension( *regAr, oddsA, *regBr, oddsB, drop_threshold );
	}

	// Do forward extension.
	closed_interval frow_interval(p.x,seqA->data.size()-1);
	dna_sequence_region_ptr regA(new dna_sequence_region(seqA, frow_interval ));
	dna_sequence_region_ptr regB(new dna_sequence_region(seqB, closed_interval(p.y, seqB->data.size() - 1) ));
	if( stop_at_repeat ) { truncate_at_repeat(*regA); truncate_at_repeat(*regB);}
	tie( max_valf, distf)  = PHMM.ungapped_forward_extension( *regA, oddsA, *regB, oddsB, drop_threshold );

	// Add extension to store.
	dstore.add( dstore.compute_diagonal(p), closed_interval(p.y-dist,p.y+distf-1) );
	return max_valf*max_val;
}

template< typename DIST >
void multi_alignment_model::set_region_dist( DIST const &rdist ) {
	dirty = true;
	typename DIST::const_iterator i = rdist.begin();
	double sum_check = 0.0;

	// For now we ignore random distribution pr to maintain score compatability with cape.
	// Besides, we may need this still to kill extension tails.
	if( i == rdist.end() ) throw std::runtime_error("multi_alignment_model: model distribution invalid");
	double r_p = *i;
	++i;

	// Now set rest. We add r_p/nmodels to each so that we sum to 1.0.
	for( int j = 0; j < nmodels; ++j, ++i ) {
		if( i == rdist.end() ) throw std::runtime_error("multi_alignment_model: model distribution invalid");
		int ms = kstates + SUBSIZE*j;
		PHMM.set_transition(state_switch,ms+state_S,*i + r_p/(double)nmodels);
		sum_check += *i + r_p/(double)nmodels;
	}
	if( fabs( 1.0 - sum_check) > 0.001 ) throw std::runtime_error("multi_alignment_model: model distribution does not sum to one");
}

#endif
