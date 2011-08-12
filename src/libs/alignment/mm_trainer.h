#ifndef MM_TRAINER_H__
#define MM_TRAINER_H__

#include "ereal.h"
#include "pair_hmm.h"
#include "dna_alignment_sequence.h"
#include "dna_sequence_region.h"
#include "nw_model_parameters.h"
#include "zero_order_background.h"
#include "fragment.h"
#include <btl/logspace.h>
#include <iostream>
#include <list>

#include "point.h"

class mm_trainer {
protected:
	typedef unsigned int uint;

	int nmodels;

	typedef pair_hmm< dna_alpha, ereal > APHMM;
	APHMM PHMM;

	enum states { state_switch = 0, state_RIa, state_RIb, kstates };
	enum sub_stats { state_S = 0, state_M, state_Ia, state_Ib, state_E, SUBSIZE };

	bool dirty;

public:
	mm_trainer( int i );

	void start_training();

	ereal end_training();

	template<typename ITERATOR>
	ereal constrained_bw_train( ITERATOR i, ITERATOR end, int nthreads );

//	template<typename ITERATOR>
//	ereal constrained_bw_train_set( ITERATOR i, ITERATOR end, int nthreads );

//	template<typename D1, typename D2, typename OUTITER>
//	void forward_anchor( D1 A, D2 B, OUTITER oit, int min_ext = 50, int max_ext = 2000 );

//	template<typename D1, typename D2, typename OUTITER>
//	void forward_anchor_b( D1 A, D2 B, OUTITER oit, int step = 50, int min_ext = 100, int max_ext = 2000, double ext_threshold = 10.0 );

//	template<typename D1, typename D2, typename O1, typename O2, typename OUTPUT_ITER>
//	void forward_anchor_c( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, OUTPUT_ITER oit, double drop_threshold );

	int num_models() const;
	nw_model_parameters_ptr get_parameters( int k ) const;
	void set_parameters( int k, nw_model_parameters const &p, int t );
	void set_default_parameters( blitz::TinyVector<double,4> q );


	template< typename DIST >
	void set_region_dist( DIST const &rdist );

	double get_region_pr( int k ) const;
	double get_region_length( int k ) const;

	double get_random_pr() const;
	double get_random_length_a() const;
	double get_random_length_b() const;

	void set_random_param( int i, int j, blitz::TinyVector<double,4> q );
	blitz::TinyVector<double,4> get_random_distribution() const;

	void set_initial_background( blitz::TinyVector<double,4> q );
	void disable_random_model();

	void set_interface_size( int i );

	void pretty_print_parameters( std::ostream &out ) const;
	void write_parameters( std::ostream &out ) const;
	void read_parameters( std::istream &in );


	//template<typename SEQA, typename SEQB, typename DP, typename RESULT, typename ODDS, typename ODDSA, typename ODDSB>
	//inline void forward_anchor_fill_dpv1( SEQA &A, SEQB &B, DP &f, RESULT &result, ODDS &odds, ODDSA &odds_data_a, ODDSB &odds_data_b, closed_interval n, closed_interval m );



};

template<typename ITER>
ereal mm_trainer::constrained_bw_train( ITER start, ITER end, int nthreads = 1 ) {

	// compute background (we assume parameters were just set so that each insdel is set to q
	//logint pback= 1.0;
	//for( ITER j = start; j != end; ++j ) {
	//	for( int i = 0; i < j->first->data.size(); ++i ) pback *= PHMM.get_emmission(IA,dna_alpha::to_index(j->first->data[i]),APHMM::NOTHING);
	//	for( int i = 0; i < j->second->data.size(); ++i ) pback *= PHMM.get_emmission(IA,dna_alpha::to_index(j->second->data[i]),APHMM::NOTHING);
	//}
	if( dirty ) { PHMM.optimize(); dirty = false; }
	return PHMM.constrained_train_bw_mt<dna_sequence_region,dna_sequence_region>(nthreads, start,end,40,true);
}

/*
template<typename D1, typename D2, typename O1, typename O2, typename OUTPUT_ITER>
void mm_trainer::forward_anchor_c( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, OUTPUT_ITER oit, double drop_threshold ) {
	using namespace std;

	// Get a list of extension points.
	list<point> extension_points;
	point last_point = PHMM.forward_extension( A, oddsA, B, oddsB, drop_threshold, back_inserter(extension_points), 100 );

	// Build a set of regions from out anchor points and output to oit.
	closed_interval intA(0,-1), intB(0,-1);
	for( list<point>::const_iterator i = extension_points.begin(); i != extension_points.end(); ++i ) {

		// break if into tail
		//if( intA.a > last_point.x || intB.a > last_point.y) break;

		//if( i->x < last_point.x && i->y < last_point.y ) {
			intA.b = i->x;
			intB.b = i->y;
		//} else {
		//	intA.b = last_point.x;
		//	intB.b = last_point.y;
		//}

		// constrain to tail=

		fragment_ptr f( new fragment() );
		f->add( intA ); f->add( intB );
		*oit = f;
		++oit;

		// Setup left bounds for next interval.
		intA.a = intA.b + 1;
		intB.a = intB.b + 1;
	}

}
*/
/*
template<typename D1, typename D2, typename OUTPUT_ITER>
void mm_trainer::forward_anchor_b( D1 A, D2 B, OUTPUT_ITER oit, int step, int min_ext, int max_ext, double ext_threshold ) {
	using namespace std;
		using namespace blitz;
		int n = A.data.size()-1, m = B.data.size()-1;

		//double drop_threshold = PHMM.get_transition(state_M,state_Ia).as_base() + 100.0*PHMM.get_transition(state_Ia,state_Ia).as_base();
		double drop_threshold = -70.0;
		//PHMM.optimize();

		// declare the data here so we don't have to keep reallocating it
		APHMM::data f(max_ext+1,max_ext+1,PHMM.M);

		// adjust initial range
		//closed_interval srange_a(0, min(max_ext-1,n) );
		//closed_interval srange_b(0, min(max_ext-1,m) );

		int total_a = 0, total_b = 0;

		// disable this for now, use a score drop like threshold instead
		min_ext = max_ext - 1;
		closed_interval srange_a(0, min(min_ext,n) );
		closed_interval srange_b(0, min(min_ext,m) );

		APHMM::fa_result result;
		bool use_fast = true;
		while(1) {
			D1 subA = A.sub_region(srange_a);
			D2 subB = B.sub_region(srange_b);

			clock_t start_time = clock();
			result = PHMM.compute_forward_anchor_b( f, odds, subA.data, subB.data, step, drop_threshold, ext_threshold, use_fast);
			clock_t end_time = clock();

			// compute background
			cerr << subA.data.range << "\t" << subB.data.range << "\t";
			cerr << result.score.as_base() << "\t" << result.max_i << "\t" << result.max_j << "\t";
			cerr << result.best_score.as_base() << "\t" << result.best_i << "\t" << result.best_j << endl;

			bool poor_extension = result.score.as_base() < ext_threshold || result.max_i <= 1 || result.max_j <= 1;

			//if( poor_extension && use_fast ) {
				//cerr << "slice_vt_fail\t" << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << "\t" << max(result.max_i,result.max_j) << endl;
				//cerr << "switch to slower algorithm" << endl;
				//use_fast = false;
			//} else if( poor_extension ) {
			if( poor_extension ) {
				cerr << "slice_vt_fail\t" << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << "\t" << max(result.max_i,result.max_j) << endl;

				// failed, just end it
				srange_a.b = srange_a.a + result.best_i-1;
				srange_b.b = srange_b.a + result.best_j-1;
				fragment_ptr f( new fragment() );
				f->add( srange_a ); f->add( srange_b );
				*oit = f;
				break;

			} else {
				if( use_fast )
				cerr << "slice_vt_good\t" << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << "\t" << max(result.max_i,result.max_j) << endl;
				else
				cerr << "slice_fw_good\t " << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << "\t" << max(result.max_i,result.max_j) << endl;

				// good block, move forward
				srange_a.b = srange_a.a + result.max_i-2;
				srange_b.b = srange_b.a + result.max_j-2;
				fragment_ptr f( new fragment() );
				f->add( srange_a ); f->add( srange_b );
				*oit = f;
				++oit;

				total_a += srange_a.size();
				total_b += srange_b.size();
				min_ext = max( max(total_a,total_b), min_ext);

				srange_a.a += result.max_i-1;
				srange_a.b = min( srange_a.a+min(max_ext-1,min_ext), n );

				srange_b.a += result.max_j-1;
				srange_b.b = min( srange_b.a+min(max_ext-1,min_ext), m );
			}
		}

}

template<typename D1, typename D2, typename OUTPUT_ITER>
void mm_trainer::forward_anchor( D1 A, D2 B, OUTPUT_ITER oit, int min_ext, int max_ext ) {
	using namespace std;
	using namespace blitz;
	int n = A.data.size()-1, m = B.data.size()-1;


	PHMM.optimize();

	// adjust initial range
	int ext = min_ext;
	closed_interval srange_a(0, min(ext-1,n) );
	closed_interval srange_b(0, min(ext-1,m) );

	APHMM::border left, top;
	APHMM::fa_result result;
	zero_order_background odds;
	odds.set_q(get_random_distribution());

	while(1) {
		D1 subA = A.sub_region(srange_a);
		D2 subB = B.sub_region(srange_b);
		cerr << subA.data.range << "\t" << subB.data.range << "\t";

		result = PHMM.compute_forward_anchor(odds, subA.data, subB.data );

		// compute background
		cerr << result.score.as_base() << "\t" << result.max_i << "\t" << result.max_j << "\t";
		cerr << result.best_score.as_base() << "\t" << result.best_i << "\t" << result.best_j << endl;

		closed_interval store_a = A.data.range, store_b = B.data.range;

		bool poor_extension = result.score.as_base() < 0 || result.max_i <= 1 || result.max_j <= 1;
		if( poor_extension && ext*2 <= max_ext ) {
			// try to increase ext
			ext *= 2;
			srange_a.b = min( srange_a.a + ext-1, n );
			srange_b.b = min( srange_b.a + ext-1, m );

		} else if( poor_extension ) {
			// failed, just end it
			srange_a.b = srange_a.a + result.best_i-1;
			srange_b.b = srange_b.a + result.best_j-1;
			fragment_ptr f( new fragment() );
			f->add( srange_a ); f->add( srange_b );
			*oit = f;
			break;

		} else {
			// good block, reset ext and move forward
			ext = min_ext;
			srange_a.b = srange_a.a + result.max_i-2;
			srange_b.b = srange_b.a + result.max_j-2;
			fragment_ptr f( new fragment() );
			f->add( srange_a ); f->add( srange_b );
			*oit = f;
			++oit;

			srange_a.a += result.max_i-1;
			srange_a.b = min( srange_a.a+ext-1, n );

			srange_b.a += result.max_j-1;
			srange_b.b = min( srange_b.a+ext-1, m );

			// end of sequence!
			if( A.data.range.size() == 0 || B.data.range.size() == 0 ) break;
		}
	}
}
*/

template< typename DIST >
void mm_trainer::set_region_dist( DIST const &rdist ) {
	dirty = true;
	typename DIST::const_iterator i = rdist.begin();

	double sum_check = 0.0;

	// first set random distribution
	if( i == rdist.end() ) throw std::runtime_error("mm_trainer: model distribution invalid");
	PHMM.set_transition(state_switch,state_RIa,*i);

	sum_check += *i;
	++i;

	// now set rest
	for( int j = 0; j < nmodels; ++j, ++i ) {
		if( i == rdist.end() ) throw std::runtime_error("mm_trainer: model distribution invalid");
		int ms = kstates + SUBSIZE*j;
		PHMM.set_transition(state_switch,ms+state_S,*i);
		sum_check += *i;
	}

	if( fabs( 1.0 - sum_check) > 0.001 ) throw std::runtime_error("mm_trainer: model distribution does not sum to one");
}

/*
template<typename SEQA, typename SEQB, typename DP, typename RESULT, typename ODDS, typename ODDSA, typename ODDSB>
inline void mm_trainer::forward_anchor_fill_dpv1( SEQA &A, SEQB &B, DP &f, RESULT &result, ODDS &odds, ODDSA &odds_data_a, ODDSB &odds_data_b, closed_interval n, closed_interval m ) {
	using namespace std;
	using namespace blitz;
	using namespace boost;
	Range all = Range::all();

	// set everything in this block to zero
	f(Range(n.a,n.b),Range(m.a,m.b),all) = ereal(0.0);

	ereal score, t;
	for( int j = m.a; j <= m.b; ++j ) {
		for( int i = n.a; i <= n.b; ++i ) {
			int x = (int)A[i-1], y = (int)B[j-1];
			t = f(i-1,j-1,get<2>(op_list(x,y)[0]))*get<3>(PHMM.op_list(x,y)[0]);
			if( f(i,j,get<2>(op_list(x,y)[0])) < t ) f(i,j,get<2>(op_list(x,y)[0])) = t;
			t = f(i-1,j-1,get<2>(op_list(x,y)[1]))*get<3>(PHMM.op_list(x,y)[1]);
			if( f(i,j,get<2>(op_list(x,y)[1])) < t ) f(i,j,get<2>(op_list(x,y)[1])) = t;
			t = f(i-1,j-1,get<2>(op_list(x,y)[2]))*get<3>(PHMM.op_list(x,y)[2]);
			if( f(i,j,get<2>(op_list(x,y)[2])) < t ) f(i,j,get<2>(op_list(x,y)[2])) = t;
			t = f(i-1,j-1,get<2>(op_list(x,y)[3]))*get<3>(PHMM.op_list(x,y)[3]);
			if( f(i,j,get<2>(op_list(x,y)[3])) < t ) f(i,j,get<2>(op_list(x,y)[3])) = t;
			t = f(i-1,j,get<2>(op_list(x,y)[4]))*get<3>(PHMM.op_list(x,y)[4]);
			if( f(i,j,get<2>(op_list(x,y)[4])) < t ) f(i,j,get<2>(op_list(x,y)[4])) = t;
			t = f(i-1,j,get<2>(op_list(x,y)[5]))*get<3>(PHMM.op_list(x,y)[5]);
			if( f(i,j,get<2>(op_list(x,y)[5])) < t ) f(i,j,get<2>(op_list(x,y)[5])) = t;
			t = f(i-1,j,get<2>(op_list(x,y)[6]))*get<3>(PHMM.op_list(x,y)[6]);
			if( f(i,j,get<2>(op_list(x,y)[6])) < t ) f(i,j,get<2>(op_list(x,y)[6])) = t;
			t = f(i,j-1,get<2>(op_list(x,y)[7]))*get<3>(PHMM.op_list(x,y)[7]);
			if( f(i,j,get<2>(op_list(x,y)[7])) < t ) f(i,j,get<2>(op_list(x,y)[7])) = t;
			t = f(i,j-1,get<2>(op_list(x,y)[8]))*get<3>(PHMM.op_list(x,y)[8]);
			if( f(i,j,get<2>(op_list(x,y)[8])) < t ) f(i,j,get<2>(op_list(x,y)[8])) = t;
			t = f(i,j-1,get<2>(op_list(x,y)[9]))*get<3>(PHMM.op_list(x,y)[9]);
			if( f(i,j,get<2>(op_list(x,y)[9])) < t ) f(i,j,get<2>(op_list(x,y)[9])) = t;
			t = f(i,j,get<2>(op_list(x,y)[10]))*get<3>(PHMM.op_list(x,y)[10]);
			if( f(i,j,get<2>(op_list(x,y)[10])) < t ) f(i,j,get<2>(op_list(x,y)[10])) = t;
			t = f(i,j,get<2>(op_list(x,y)[11]))*get<3>(PHMM.op_list(x,y)[11]);
			if( f(i,j,get<2>(op_list(x,y)[11])) < t ) f(i,j,get<2>(op_list(x,y)[11])) = t;
			t = f(i,j,get<2>(op_list(x,y)[12]))*get<3>(PHMM.op_list(x,y)[12]);
			if( f(i,j,get<2>(op_list(x,y)[12])) < t ) f(i,j,get<2>(op_list(x,y)[12])) = t;
			t = f(i,j,get<2>(op_list(x,y)[13]))*get<3>(PHMM.op_list(x,y)[13]);
			if( f(i,j,get<2>(op_list(x,y)[13])) < t ) f(i,j,get<2>(op_list(x,y)[13])) = t;
			t = f(i,j,get<2>(op_list(x,y)[14]))*get<3>(PHMM.op_list(x,y)[14]);
			if( f(i,j,get<2>(op_list(x,y)[14])) < t ) f(i,j,get<2>(op_list(x,y)[14])) = t;


			score = f(i,j,end_state)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
			if( score > result.best_score ) {
				result.best_i = i;
				result.best_j = j;
				result.best_score = score;
			}
		}
	}
}
*/
#endif
