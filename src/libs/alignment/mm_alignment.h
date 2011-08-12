#ifndef MM_ALIGNMENT_H__
#define MM_ALIGNMENT_H__
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <btl/logspace.h>
#include "pairwise_dna_alignment.h"
#include "dna_sequence_region.h"
#include "nw_model_parameters.h"
#include "annotation.h"
#include "ereal.h"


class mm_alignment {
protected:
	typedef unsigned int uint;

	struct dpdata {
		ereal s;
		uint16_t from;
		dpdata() : s( 0.0 ), from( state_invalid ) {}
		void reset() { s = 0.0; from = state_invalid; }
	};
	enum states { state_invalid = 0, state_RIa, state_RIb, state_switch, state_S, state_M, state_Ia, state_Ib, state_E };
	enum region_types { RANDOM_MODEL = 0, USER_MODEL };

	struct dpcell {
		dpdata S;		// start state
		dpdata E;		// end state
		dpdata M;		// match state
		dpdata Ia;		// insert into A
		dpdata Ib;		// insert into B
		void reset() { S.reset(); E.reset(); M.reset(); Ia.reset(); Ib.reset(); }
	};

	struct dprandomcell {
		dpdata Ia;		// random insert A
		dpdata Ib;		// random insert b
		void reset() { Ia.reset(); Ib.reset(); }
	};

	typedef blitz::Array<dpcell,3> dpmodel;
	typedef blitz::Array<dpdata,2> dpjoin;
	typedef blitz::Array<dprandomcell,2> dprandom;


	struct random_model_parameters {
		blitz::TinyVector<ereal,4>	q;
		ereal p_a;
		ereal p_1ma;
		ereal p_b;
		ereal p_1mb;
		random_model_parameters();
		void set( int i, int j, blitz::TinyVector<double,4> q );
	};

	struct model_parameters {
		blitz::Array<ereal,2>			p;	// substitution matrix
		blitz::TinyVector<ereal,4>	q;	// background
		ereal p_at;							// transition probabilities
		ereal p_bt;
		ereal p_1mbt;
		ereal p_1m2at;
		ereal p_a;
		ereal p_1m2a;
		ereal p_1mt;
		nw_model_parameters nwmp;
		void set( nw_model_parameters const &m, double t );
		model_parameters();
		model_parameters( model_parameters const &other );
	};
	typedef std::vector<model_parameters> model_parameter_set;
	typedef std::vector<ereal > model_distribution;

	struct model_data {
		// dp matrix for all models
		dpmodel m;

		// dp matrix for switch state
		dpjoin  s;

		dprandom r;
	} dp, fw, t, e;

	// parameters for the random model
	random_model_parameters rp;

	// parameters for each model
	model_parameter_set mp;

	// Parameters for choosing a model. This must sum to 1 - p_rand.
	model_distribution p_sw;

	// probability of choosing the random model
	ereal p_rand;

	// number of models
	int nmodels;

	// interface size for constrained dp
	int interface;

	/** Do we need to rescale scores during alignments. */
	bool rescale_scores;

	// data for interface
	struct alignment_chain {
		int sid;										// in a constrained alignment, this is the segment the chain came from
		pairwise_dna_alignment align;					// alignment
		annotation_ptr ann;								// annotation

		boost::shared_ptr<alignment_chain> prev;		// previous chain
		int prev_i;										// index of previous chain, zero if no previous chain
		int prev_k;										// previous model
		int prev_state;									// start in model
		int i;											// start location in previous chain
		int j;
		int i_end;										// end of chain index
		int j_end;
	};
	typedef boost::shared_ptr<alignment_chain> alignment_chain_ptr;

	struct interface_data {
		// store model data
		model_data dp;

		// chain for each state in each model
		std::vector< alignment_chain_ptr > alignments;

		interface_data( int k, int i ) {
			dp.m.resize(k,i*2,1);
			dp.s.resize(i*2,1);
			dp.r.resize(i*2,1);
			alignments.resize((k*3+2)*i*2); // make sure we have enough for (M,Ia,Ib (3) * k  + RIa + RIb (2)) * interface size alignments
		}
	};

	inline unsigned int map_to( int i, int k, int state ) {
		switch( state ) {
			case state_M:
				return (2 + 3*nmodels)*i + 2 + 3*k;
				break;
			case state_Ia:
				return (2 + 3*nmodels)*i + 2 + 3*k + 1;
				break;
			case state_Ib:
				return (2 + 3*nmodels)*i + 2 + 3*k + 2;
				break;
			case state_RIa:
				return (2 + 3*nmodels)*i;
				break;
			case state_RIb:
				return (2 + 3*nmodels)*i + 1;
				break;
		}
		std::cerr << state << std::endl;
		assert( 0 && "invalid state" );
		return 0;
	}

	// for training
	typedef blitz::Array<int,3> pair_count;
	typedef blitz::Array<int,2> base_count;

	struct all_counts {
		pair_count pairs;
		std::vector<int> gaps_a;
		std::vector<int> gaps_b;
		std::vector<int> gap_lengths_a;
		std::vector<int> gap_lengths_b;
		std::vector<int> region_uses;
		std::vector<int> region_lengths;

		base_count bases_a;
		base_count bases_b;

		all_counts() {}

		void reset( int k ) {
			// adjust for other models
			k = USER_MODEL + k;
			pairs.resize(k,dna_alpha::SIZE,dna_alpha::SIZE);
			pairs = 0;
			gaps_a.clear();
			gaps_a.resize(k,0);
			gaps_b.clear();
			gaps_b.resize(k,0);
			gap_lengths_a.clear();
			gap_lengths_a.resize(k,0);
			gap_lengths_b.clear();
			gap_lengths_b.resize(k,0);
			region_uses.clear();
			region_uses.resize(k,0);
			region_lengths.clear();
			region_lengths.resize(k,0);
			bases_a.resize(k,dna_alpha::SIZE);
			bases_a = 0;
			bases_b.resize(k,dna_alpha::SIZE);
			bases_b = 0;
		}
	} my_counts;

	double my_pr;

	pairwise_dna_alignment viterbi_decode( dna_sequence_region_data &a, dna_sequence_region_data &b );
	void viterbi_init_borders( model_data &mydp, model_data &init, dna_sequence_region_data &a, dna_sequence_region_data &b );
	void viterbi_filldp( dna_sequence_region_data &a, dna_sequence_region_data &b );
	alignment_chain_ptr viterbi_traceback( int k1, int k2, int i, int j, dna_sequence_region_data &a, dna_sequence_region_data &b, uint16_t p = 0, bool chain_required = false , int s = -1);

	ereal p( int k, int i, int j, dna_sequence_region_data &a, dna_sequence_region_data &b );
	ereal q( int k, int i, dna_sequence_region_data &a );
	ereal q_rand( int i, dna_sequence_region_data &a );

	void do_statistics( std::pair<pairwise_dna_alignment,annotation_ptr> const &r, all_counts &c );
	void update_probabilities( all_counts &counts );


	void perform_chain( alignment_chain &c, interface_data &id, interface_data &id2 );
	std::pair<pairwise_dna_alignment,annotation_ptr> collapse_chain( alignment_chain &c );


public:
	mm_alignment( int i );
	pairwise_dna_alignment align( dna_sequence_region &a, dna_sequence_region &b );
	std::pair<pairwise_dna_alignment,annotation_ptr> align_and_annotate( dna_sequence_region &a, dna_sequence_region &b );

	void start_training();
	double end_training();

	template<typename ITERATOR>
	void viterbi_train( ITERATOR i, ITERATOR end );

	template<typename ITERATOR>
	void constrained_viterbi_train( ITERATOR i, ITERATOR end );

	template<typename ITERATOR>
	std::pair<pairwise_dna_alignment,annotation_ptr> constrained_align_and_annotate( ITERATOR i, ITERATOR end );

	int num_models();
	nw_model_parameters const &get_parameters( int k );
	void set_parameters( int k, nw_model_parameters const &p, int t );

	template< typename DIST >
	void set_region_dist( DIST const &rdist );

	double get_region_pr( int k );
	double get_region_length( int k );

	double get_random_pr();
	double get_random_length_a();
	double get_random_length_b();
	void set_random_param( int i, int j, blitz::TinyVector<double,4> q );
	blitz::TinyVector<ereal,4> get_random_distribution();

	void set_default_parameters( blitz::TinyVector<double,4> q );

	void set_initial_background( blitz::TinyVector<double,4> q );
	void disable_random_model();

	void set_interface_size( int i );
	void read_parameters( std::istream &in );


};

template< typename DIST >
void mm_alignment::set_region_dist( DIST const &rdist ) {
	typename DIST::const_iterator i = rdist.begin();

	double sum_check = 0.0;

	// first set random distribution
	if( i == rdist.end() ) throw std::runtime_error("mm_alignment: model distribution invalid");
	p_rand = *i;

	sum_check += *i;
	++i;

	// now set rest
	for( int j = 0; j < nmodels; ++j, ++i ) {
		if( i == rdist.end() ) throw std::runtime_error("mm_alignment: model distribution invalid");
		p_sw[j] = *i;
		sum_check += *i;
	}

	if( fabs( 1.0 - sum_check) > 0.001 ) throw std::runtime_error("mm_alignment: model distribution does not sum to one");
}


inline ereal mm_alignment::p( int k, int i, int j, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
	return mp[k].p((int)a[i-1],(int)b[j-1]);
}

inline ereal mm_alignment::q_rand( int i, dna_sequence_region_data &a ) {
	return rp.q((int)a[i-1]);
}


inline ereal mm_alignment::q( int k, int i, dna_sequence_region_data &a ) {
	return mp[k].q((int)a[i-1]);
}

// train on a list of sequences
template<typename ITERATOR>
void mm_alignment::viterbi_train( ITERATOR it, ITERATOR end ) {
	using namespace blitz;

	all_counts counts(nmodels);

	for( ; it != end; ++it ) {
		std::pair<pairwise_dna_alignment,annotation_ptr> r = align_and_annotate( *(it->first), *(it->second) );
		do_statistics( r, counts );
	}

	update_probabilities( counts );

}

// train on a list of sequences represeting a large constrained alignment
template<typename ITERATOR>
void mm_alignment::constrained_viterbi_train( ITERATOR it, ITERATOR end ) {
	using namespace blitz;

	std::pair<pairwise_dna_alignment,annotation_ptr> r = constrained_align_and_annotate( it, end );
	do_statistics( r, my_counts );
	my_pr += r.first.score;
}



// perform a constrained alignment
template<typename ITERATOR>
std::pair<pairwise_dna_alignment,annotation_ptr> mm_alignment::constrained_align_and_annotate( ITERATOR start, ITERATOR end ) {
	ITERATOR next = start, it = start;
	++next;

	// initialize interface with starting parameters
	interface_data idata(nmodels,interface), idata_next(nmodels,interface);
	model_data null;

	// need final region size for later
	int n = 0, m = 0, interface_size = interface*2;

	std::vector< ereal > score_adjustments;

	int segment = 0;
	for( ; it != end; ++it, ++next ) {
		// first adjust region to align based on interface
		dna_sequence_region region_a = *(it->first), region_b = *(it->second);

		// if not at end, adjust interface
		if( next != end ) region_a.data.range.b += interface;

		if( it != start ) region_a.data.range.a -= (interface-1);
#ifdef DEBUG_OUTPUT
		std::cerr << "Aligning " << region_a.data.range << " vs " << region_b.data.range << std::endl;
#endif
		if( region_a.data.size() < (unsigned int)interface_size ) interface_size = (int) region_a.data.size();

		// allocate dp matricies and fill
		n = (int) region_a.data.size(), m = (int)region_b.data.size();

		// allocate dp matricies if needed
		if( dp.m.extent(blitz::secondDim) < (n+1) || dp.m.extent(blitz::thirdDim) < (m+1) ) dp.m.resize(nmodels,n+1,m+1);
		if( dp.s.extent(blitz::firstDim) < (n+1) || dp.s.extent(blitz::secondDim) < (m+1) ) dp.s.resize(n+1,m+1);
		if( dp.r.extent(blitz::firstDim) < (n+1) || dp.r.extent(blitz::secondDim) < (m+1) ) dp.r.resize(n+1,m+1);

		// zero matrix borders
		for( int i = 0; i <= n; ++i ) {
			for( int k = 0; k < nmodels; ++k ) dp.m(k,i,0).reset();
			dp.s(i,0).reset();
			dp.r(i,0).reset();
		}
		for( int j = 0; j <= m; ++j ) {
			for( int k = 0; k < nmodels; ++k ) dp.m(k,0,j).reset();
			dp.s(0,j).reset();
			dp.r(0,j).reset();
		}

		if( it == start ) viterbi_init_borders( dp, null, region_a.data, region_b.data );
		else viterbi_init_borders( dp, idata.dp, region_a.data, region_b.data );
		viterbi_filldp(region_a.data, region_b.data);

		// find all optimal alignment chains
		for( int j = 0, i = 1; i <= interface_size; ++i ) {
			for( int k = 0; k < nmodels; ++k ) {
				j = map_to( i-1, k, state_M );
				idata_next.alignments[j] = viterbi_traceback(k,state_M,n-interface_size+i,m,region_a.data,region_b.data,j+1, it != start, segment );
				perform_chain( *(idata_next.alignments[j]), idata, idata_next );
				j = map_to( i-1, k, state_Ia );
				idata_next.alignments[j] = viterbi_traceback(k,state_Ia,n-interface_size+i,m,region_a.data,region_b.data,j+1, it != start, segment );
				perform_chain( *(idata_next.alignments[j]), idata, idata_next );
				j = map_to( i-1, k, state_Ib );
				idata_next.alignments[j] = viterbi_traceback(k,state_Ib,n-interface_size+i,m,region_a.data,region_b.data,j+1, it != start, segment );
				perform_chain( *(idata_next.alignments[j]), idata, idata_next );
			}
			j = map_to( i-1, 0, state_RIa );
			idata_next.alignments[j] = viterbi_traceback(-1,state_RIa,n-interface_size+i,m,region_a.data,region_b.data,j+1, it != start, segment );
			perform_chain( *(idata_next.alignments[j]), idata, idata_next );
			j = map_to( i-1, 0, state_RIb );
			idata_next.alignments[j] = viterbi_traceback(-1,state_RIb,n-interface_size+i,m,region_a.data,region_b.data,j+1, it != start, segment );
			perform_chain( *(idata_next.alignments[j]), idata, idata_next );

		}

		// copy border data back to idata
		for( int i = 0; i < interface_size; ++i ) {
			for( int k = 0; k < nmodels; ++k ) idata_next.dp.m(k,i,0) = dp.m(k,n-interface_size+1+i,m);
			idata_next.dp.s(i,0) = dp.s(n-interface_size+1+i,m);
			idata_next.dp.r(i,0) = dp.r(n-interface_size+1+i,m);
		}

		// Scores can grow and overflow, so we adjust them downwards on each interval.
		// First we find the (largest probability) and set that to 1. This is equivalent to setting the
		// largest score to 0.
		if( rescale_scores ) {
			ereal max_pr = 0.0;
			for( int i = 0; i < interface_size; ++i ) {
				for( int k = 0; k < nmodels; ++k ) {
					if( idata_next.dp.m(k,i,0).M.s > max_pr ) max_pr = idata_next.dp.m(k,i,0).M.s;
					if( idata_next.dp.m(k,i,0).Ia.s > max_pr ) max_pr = idata_next.dp.m(k,i,0).Ia.s;
					if( idata_next.dp.m(k,i,0).Ib.s > max_pr ) max_pr = idata_next.dp.m(k,i,0).Ib.s;
					if( idata_next.dp.m(k,i,0).S.s > max_pr ) max_pr = idata_next.dp.m(k,i,0).S.s;
					if( idata_next.dp.m(k,i,0).E.s > max_pr ) max_pr = idata_next.dp.m(k,i,0).E.s;
				}
				if( idata_next.dp.s(i,0).s > max_pr ) max_pr = idata_next.dp.s(i,0).s;
				if( idata_next.dp.r(i,0).Ia.s > max_pr ) max_pr = idata_next.dp.r(i,0).Ia.s;
				if( idata_next.dp.r(i,0).Ib.s > max_pr ) max_pr = idata_next.dp.r(i,0).Ib.s;
			}

			// now divide everything by max_pr so that the max_pr becomes 1.0
			for( int i = 0; i < interface_size; ++i ) {
				for( int k = 0; k < nmodels; ++k ) {
					idata_next.dp.m(k,i,0).M.s /= max_pr;
					idata_next.dp.m(k,i,0).Ia.s /= max_pr;
					idata_next.dp.m(k,i,0).Ib.s /= max_pr;
				}
				idata_next.dp.s(i,0).s /= max_pr;
				idata_next.dp.r(i,0).Ia.s /= max_pr;
				idata_next.dp.r(i,0).Ib.s /= max_pr;
			}
			score_adjustments.push_back( max_pr );
		}
		idata = idata_next;
		++segment;
	}

	// now we need to collapse the optimal chain into one alignment
	// alignment sequences
	std::pair<pairwise_dna_alignment,annotation_ptr> r;
	if( dp.s(n,m).from == state_RIb ) {
		r = collapse_chain( *(idata.alignments[map_to(interface_size-1,-1,state_RIb)]) );
	} else {
		int k = dp.s(n,m).from - state_E;
		int state = dp.m(k,n,m).E.from;
		r = collapse_chain( *(idata.alignments[map_to(interface_size-1,k,state)]) );
	}

	// compute odds
	ereal odds = 1.0;
	for( annotation_iter i = r.second->begin(); i != r.second->end(); ++i ) {
		for( int c = i->a; c <= i->b; ++c ) {
			if( r.first.a->data[c] != dna_alignment_alpha::GAP ) odds *= rp.q((int)dna_alpha::symbol(r.first.a->data[c]));
			if( r.first.b->data[c] != dna_alignment_alpha::GAP ) odds *= rp.q((int)dna_alpha::symbol(r.first.b->data[c]));
		}
	}

	ereal final_score = dp.s(n,m).s/odds;
	for( std::vector< ereal >::iterator i = score_adjustments.begin(); i != score_adjustments.end(); ++i )
		final_score *= *i;

	r.first.score = final_score.as_base();
	return r;
}


#endif

