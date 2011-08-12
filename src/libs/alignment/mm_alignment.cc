#include "mm_alignment.h"
#include "io.h"
#include "ereal.h"

#include <boost/tuple/tuple.hpp>

mm_alignment::mm_alignment( int i ) : nmodels(i), interface(1), rescale_scores(false) {
	mp.resize(i);
	p_sw.resize(i);
	set_default_parameters(blitz::TinyVector<double,4>(0.25,0.25,0.25,0.25));
}

void mm_alignment::set_initial_background( blitz::TinyVector<double,4> q ) {
	set_default_parameters(q);
	rp.set( 500, 500, q );
}

mm_alignment::model_parameters::model_parameters() {
	p.resize(dna_alpha::SIZE,dna_alpha::SIZE);
}

mm_alignment::random_model_parameters::random_model_parameters() {
	set( 500, 500,blitz::TinyVector<double,4>(0.25,0.25,0.25,0.25) );
}

void mm_alignment::random_model_parameters::set( int a, int b, blitz::TinyVector<double,4> nq ) {
	p_a = 1.0/(1.0 + 1.0/(double)a);
	p_1ma = 1.0/(1.0+(double)a);
	p_b = 1.0/(1.0 + 1.0/(double)b);
	p_1mb = 1.0/(1.0 + (double)b);
	for( int i = 0; i < 4; ++i) q(i) = nq(i);
}


mm_alignment::model_parameters::model_parameters( mm_alignment::model_parameters const &other) {
	p.resize(dna_alpha::SIZE,dna_alpha::SIZE);
	p = other.p;
	q = other.q;
	p_at = other.p_at;
	p_1m2at = other.p_1m2at;
	p_a = other.p_a;
	p_1m2at = other.p_1m2at;
	p_bt = other.p_bt;
	p_1mbt = other.p_1mbt;
	p_1mt = other.p_1mt;
	nwmp = other.nwmp;
}


void mm_alignment::model_parameters::set( nw_model_parameters const &nwp, double t ) {
	// set open
   p_at = nwp.pr_open*t;
	p_1m2at = (1.0-2.0*nwp.pr_open)*t;
	p_a = nwp.pr_open;
	p_1m2a = 1.0 - 2.0*nwp.pr_open;

	// Since a gap is length 1 by the virtue of opening it, we adjust n down by 1.
   // The remaining gap extension probability corresponds to the mean of a
   // geometric distribution allowing zero extensions.
   double pr_b = 1.0-1.0/nwp.mean_gap_length;
   p_bt = pr_b*t;
   p_1mbt = (1.0 - pr_b)*t;

	p_1mt = 1.0 - t;

	// set background
	for( int i = 0; i < 4; i++ ) {
      q(i) = nwp.q(i);
   }


	// set p
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			p(i,j) = nwp.p(i,j);
		}
   }

	// store nw_model_parameters for information purposes
	nwmp = nwp;
}

int mm_alignment::num_models() {
	return nmodels;
}

nw_model_parameters const &mm_alignment::get_parameters( int k ) {
	return mp[k].nwmp;
}

void mm_alignment::set_parameters( int k, nw_model_parameters const &p, int t ) {
	mp[k].set(p,1.0-1.0/(double)t);
}

// disable the random model
void mm_alignment::disable_random_model() {
	// we redistribution the probability assigned to the random model to the remaining
	// models evenly
	ereal prslice = (double)p_rand/(double)nmodels;
	for( int k = 0; k < nmodels; ++k ) p_sw[k] += prslice;
	p_rand = 0.0;
}


double mm_alignment::get_random_pr() {
	return (double)p_rand;
}

double mm_alignment::get_random_length_a() {
	return 1.0/(double)rp.p_1ma;
}

double mm_alignment::get_random_length_b() {
	return 1.0/(double)rp.p_1mb;
}

void mm_alignment::set_random_param( int i, int j, blitz::TinyVector<double,4> q ) {
	rp.set(i,j,q);
}

blitz::TinyVector<ereal,4> mm_alignment::get_random_distribution() {
	return rp.q;
}

double mm_alignment::get_region_pr( int k ) {
	return (double)p_sw[k];
}

double mm_alignment::get_region_length( int k ) {
	return 1.0/(double)mp[k].p_1mt;
}

void mm_alignment::set_default_parameters( blitz::TinyVector<double,4> q ) {
	using namespace std;
	using namespace boost;
	assert(nmodels < 30);

	typedef boost::tuple<double,double,double> double3;
	vector<double3> custom_pts;
	custom_pts.push_back( double3(0.4,0.0625,2.0) );
	custom_pts.push_back( double3(0.2,0.03125,1.333) );
	custom_pts.push_back( double3(0.6,0.0625,3.41424) );
	custom_pts.push_back( double3(0.8,0.03125,2.0) );

	/*
	//static unsigned int const greedy_opt[] = {106, 88, 150, 62, 59, 77, 103, 145, 44, 121, 132, 31, 84, 93, 72,
	//131, 146, 75, 28, 54, 61, 78, 134, 138, 142, 46, 58, 63, 115, 119 };
	static unsigned int const greedy_opt[] = {16, 85, 40, 75, 59, 77, 103, 145, 44, 121, 132, 31, 84, 93, 72,
	 131, 146, 75, 28, 54, 61, 78, 134, 138, 142, 46, 58, 63, 115, 119 };

   vector<double> gap_lengths, gap_opens, distances;
   gap_opens.push_back( 0.03125 ); // -5
   gap_opens.push_back( 0.0625  ); // -4
   gap_opens.push_back( 0.125   ); // -3
   gap_opens.push_back( 0.25    ); // -2

   gap_lengths.push_back( 1.333    ); // -2
   gap_lengths.push_back( 2.0      ); // -1
   gap_lengths.push_back( 3.414214 ); // -0.5
   gap_lengths.push_back( 6.285214 ); // -0.25

	distances.push_back( 0.1 );
   distances.push_back( 0.2 );
   distances.push_back( 0.3 );
   distances.push_back( 0.4 );
   distances.push_back( 0.5 );
   distances.push_back( 0.6 );
   distances.push_back( 0.7 );
   distances.push_back( 0.8 );
   distances.push_back( 0.9 );
   distances.push_back( 1.0 );

	std::vector<nw_model_parameters> all;
	for( vector<double>::const_iterator d = distances.begin(); d != distances.end(); ++d ) {
		for( vector<double>::const_iterator a = gap_opens.begin(); a != gap_opens.end(); ++a ) {
			for( vector<double>::const_iterator L = gap_lengths.begin(); L != gap_lengths.end(); ++L ) {
				all.push_back( nw_model_parameters() );
				all.back().set_p_hky( *d, 1.0 );
				all.back().pr_open = *a;
				all.back().mean_gap_length = *L;

			}
		}
   }
	*/

	double pr_rand =  1.0/10000.0;
	p_rand = pr_rand;
	rp.set( 500, 500,q );


	assert( (unsigned int)nmodels <= custom_pts.size() );
	for( int k = 0; k < nmodels; ++k ) {
		nw_model_parameters newp;
		newp.q = q;
		newp.set_p_hky( custom_pts[k].get<0>(), 1.0 );
		newp.pr_open = custom_pts[k].get<1>();
		newp.mean_gap_length = custom_pts[k].get<2>();
		mp[k].set( newp, 1.0 - 1.0/200 );
		p_sw[k] = (1.0-pr_rand)/nmodels;
	}
}

std::pair<pairwise_dna_alignment,annotation_ptr> mm_alignment::align_and_annotate( dna_sequence_region &a, dna_sequence_region &b ) {
	int n = (int)a.data.size(), m = (int)b.data.size();

	// allocate dp matricies if needed
	if( dp.m.extent(blitz::secondDim) < (n+1) || dp.m.extent(blitz::thirdDim) < (m+1) ) dp.m.resize(nmodels,n+1,m+1);
	if( dp.s.extent(blitz::firstDim) < (n+1) || dp.s.extent(blitz::secondDim) < (m+1) ) dp.s.resize(n+1,m+1);
	if( dp.r.extent(blitz::firstDim) < (n+1) || dp.r.extent(blitz::secondDim) < (m+1) ) dp.r.resize(n+1,m+1);

	model_data null;
	viterbi_init_borders( dp, null, a.data, b.data );
	viterbi_filldp( a.data, b.data );
	alignment_chain_ptr r = viterbi_traceback(-1,state_switch,n,m,a.data,b.data);

	// compute odds
	ereal odds = 1.0;
	for( annotation_iter i = r->ann->begin(); i != r->ann->end(); ++i ) {
		for( int c = i->a; c <= i->b; ++c ) {
			if( r->align.a->data[c] != dna_alignment_alpha::GAP ) odds *= mp[i->t].q((int)dna_alpha::symbol(r->align.a->data[c]));
			if( r->align.b->data[c] != dna_alignment_alpha::GAP ) odds *= mp[i->t].q((int)dna_alpha::symbol(r->align.b->data[c]));
		}
	}
	r->align.score = (dp.s(n,m).s/odds).as_base();

	return make_pair(r->align,r->ann);
}

pairwise_dna_alignment mm_alignment::align( dna_sequence_region &a, dna_sequence_region &b ) {
	std::pair<pairwise_dna_alignment,annotation_ptr> r = align_and_annotate(a,b);
	return r.first;
}


void mm_alignment::viterbi_init_borders( model_data &mydp, model_data &init, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
	int n = (int)a.size(), m = (int)b.size();

	// init size will always be a multiple of two
	int init_size = init.m.extent(blitz::secondDim);
	assert( init_size % 2 == 0 );

	if( init_size > 0 ) {
		for( int i = 0; i < std::min(init_size,n); ++i ) {
			mydp.s(i,0) = init.s(i,0);
			mydp.r(i,0) = init.r(i,0);
		}
	} else {
		mydp.s(0,0).s = 1.0;
		mydp.r(1,0).Ia.s = mydp.s(0,0).s*p_rand*q_rand(1,a);
		mydp.r(1,0).Ia.from = state_switch;
	}

	// set random model init
	for( int i = (init_size > 0) ? init_size : 2; i <= n; ++i ) {
		mydp.r(i,0).Ia.s = mydp.r(i-1,0).Ia.s*rp.p_a*q_rand(i,a);
		mydp.r(i,0).Ia.from = state_RIa;
	}

	for( int k = 0; k < nmodels; ++k ) {
		if( init_size > 0 )  {
			for( int i = 0; i < std::min(init_size,n); ++i ) mydp.m(k,i,0) = init.m(k,i,0);
		} else {
			mydp.m(k,0,0).S.s = mydp.s(0,0).s*p_sw[k];
			mydp.m(k,0,0).S.from = state_switch;
			mydp.m(k,1,0).Ia.s = mydp.m(k,0,0).S.s*mp[k].p_a*q(k,1,a);
			mydp.m(k,1,0).Ia.from = state_S;
			mydp.m(k,1,0).E.s = mydp.m(k,1,0).Ia.s*mp[k].p_1mt;
			mydp.m(k,1,0).E.from = state_Ia;
			if( mydp.s(1,0).s < mydp.m(k,1,0).E.s ) {
				mydp.s(1,0).s = mydp.m(k,1,0).E.s;
				mydp.s(1,0).from = state_E + k;
			}
		}

		for( int i = (init_size > 0) ? init_size : 2; i <= n; ++i ) {
			mydp.m(k,i,0).Ia.s = mydp.m(k,i-1,0).Ia.s*mp[k].p_bt*q(k,i,a);
			mydp.m(k,i,0).Ia.from = state_Ia;
			mydp.m(k,i,0).E.s = mydp.m(k,i,0).Ia.s*mp[k].p_1mt;
			mydp.m(k,i,0).E.from = state_Ia;
			if( mydp.s(i,0).s < mydp.m(k,i,0).E.s ) {
				mydp.s(i,0).s = mydp.m(k,i,0).E.s;
				mydp.s(i,0).from = state_E + k;
			}
		}

		if( m > 0 ) {
			mydp.m(k,0,1).Ib.s = mydp.m(k,0,0).S.s*mp[k].p_a*q(k,1,b);
			mydp.m(k,0,1).Ib.from = state_S;
			// if we have initialization data we must consider coming from an open gap
			if( init_size > 0 ) {
				ereal sc = mydp.m(k,0,0).Ib.s*mp[k].p_bt*q(k,1,b);
				if( sc > mydp.m(k,0,1).Ib.s ) {
					mydp.m(k,0,1).Ib.s = sc;
					mydp.m(k,0,1).Ib.from = state_Ib;
				}
			}
			mydp.m(k,0,1).E.s = mydp.m(k,0,1).Ib.s*mp[k].p_1mt;
			mydp.m(k,0,1).E.from = state_Ib;
			if( mydp.s(0,1).s < mydp.m(k,0,1).E.s ) {
				mydp.s(0,1).s = mydp.m(k,0,1).E.s;
				mydp.s(0,1).from = state_E + k;
			}
		}

		for( int j = 2; j <= m; ++j ) {
			mydp.m(k,0,j).Ib.s = mydp.m(k,0,j-1).Ib.s*mp[k].p_bt*q(k,j,b);
			mydp.m(k,0,j).Ib.from = state_Ib;
			mydp.m(k,0,j).E.s = mydp.m(k,0,j-1).Ib.s*mp[k].p_1mt;
			mydp.m(k,0,j).E.from = state_Ib;
			if( mydp.s(0,j).s < mydp.m(k,0,j).E.s ) {
				mydp.s(0,j).s = mydp.m(k,0,j).E.s;
				mydp.s(0,j).from = state_E + k;
			}
		}
	}
}



void mm_alignment::viterbi_filldp( dna_sequence_region_data &a, dna_sequence_region_data &b ) {
	int n = (int)a.size(), m = (int)b.size();
	using namespace std;

	ereal sc;
	for( int i = 1; i <= n; ++i ) {
	for( int j = 1; j <= m; ++j ) {
		dp.s(i,j).s = 0.0;

		// do random model
		dp.r(i,j).Ia.s = dp.r(i-1,j).Ia.s*rp.p_a*q_rand(i,a);
		dp.r(i,j).Ia.from = state_RIa;
		sc = dp.s(i-1,j).s*p_rand*q_rand(i,a);
		if( sc > dp.r(i,j).Ia.s ) { dp.r(i,j).Ia.s = sc; dp.r(i,j).Ia.from = state_switch; }

		dp.r(i,j).Ib.s = dp.r(i,j-1).Ib.s*rp.p_b*q_rand(j,b);
		dp.r(i,j).Ib.from = state_RIb;
		sc = dp.r(i,j-1).Ia.s*rp.p_1ma*q_rand(j,b);
		if( sc > dp.r(i,j).Ib.s ) { dp.r(i,j).Ib.s = sc; dp.r(i,j).Ib.from = state_RIa; }

		dp.s(i,j).s = dp.r(i,j).Ib.s*rp.p_1mb;
		dp.s(i,j).from = state_RIb;

		for( int k = 0; k < nmodels; ++k ) {
			// state M
			dp.m(k,i,j).M.s = dp.m(k,i-1,j-1).S.s*mp[k].p_1m2at*p(k,i,j,a,b);
			dp.m(k,i,j).M.from = state_S;

			sc = dp.m(k,i-1,j-1).M.s*mp[k].p_1m2at*p(k,i,j,a,b);
			if( sc > dp.m(k,i,j).M.s ) { dp.m(k,i,j).M.s = sc; dp.m(k,i,j).M.from = state_M; }

			sc = dp.m(k,i-1,j-1).Ia.s*mp[k].p_1mbt*p(k,i,j,a,b);
			if( sc > dp.m(k,i,j).M.s ) { dp.m(k,i,j).M.s = sc; dp.m(k,i,j).M.from = state_Ia; }

			sc = dp.m(k,i-1,j-1).Ib.s*mp[k].p_1mbt*p(k,i,j,a,b);
			if( sc > dp.m(k,i,j).M.s ) { dp.m(k,i,j).M.s = sc; dp.m(k,i,j).M.from = state_Ib; }

			// state Ia
			dp.m(k,i,j).Ia.s = dp.m(k,i-1,j).S.s*mp[k].p_a*q(k,i,a);
			dp.m(k,i,j).Ia.from = state_S;

			sc = dp.m(k,i-1,j).M.s*mp[k].p_at*q(k,i,a);
			if( sc > dp.m(k,i,j).Ia.s ) { dp.m(k,i,j).Ia.s = sc; dp.m(k,i,j).Ia.from = state_M; }

			sc = dp.m(k,i-1,j).Ia.s*mp[k].p_bt*q(k,i,a);
			if( sc > dp.m(k,i,j).Ia.s ) { dp.m(k,i,j).Ia.s = sc; dp.m(k,i,j).Ia.from = state_Ia; }

			// state Ib
			dp.m(k,i,j).Ib.s = dp.m(k,i,j-1).S.s*mp[k].p_a*q(k,j,b);
			dp.m(k,i,j).Ib.from = state_S;

			sc = dp.m(k,i,j-1).M.s*mp[k].p_at*q(k,j,b);
			if( sc > dp.m(k,i,j).Ib.s ) { dp.m(k,i,j).Ib.s = sc; dp.m(k,i,j).Ib.from = state_M; }

			sc = dp.m(k,i,j-1).Ib.s*mp[k].p_bt*q(k,j,b);
			if( sc > dp.m(k,i,j).Ib.s ) { dp.m(k,i,j).Ib.s = sc; dp.m(k,i,j).Ib.from = state_Ib; }

			// end state
			dp.m(k,i,j).E.s = dp.m(k,i,j).Ia.s*mp[k].p_1mt;
			dp.m(k,i,j).E.from = state_Ia;

			sc = dp.m(k,i,j).M.s*mp[k].p_1mt;
			if( sc > dp.m(k,i,j).E.s ) { dp.m(k,i,j).E.s = sc; dp.m(k,i,j).E.from = state_M; }

			sc = dp.m(k,i,j).Ib.s*mp[k].p_1mt;
			if( sc > dp.m(k,i,j).E.s ) { dp.m(k,i,j).E.s = sc; dp.m(k,i,j).E.from = state_Ib; }

			// switch state
			if( dp.m(k,i,j).E.s > dp.s(i,j).s ) {
				dp.s(i,j).s = dp.m(k,i,j).E.s;
				dp.s(i,j).from = state_E + k;
			}
		} // k

		for( int k = 0; k < nmodels; ++k ) {
			// start state
			dp.m(k,i,j).S.s = dp.s(i,j).s*p_sw[k];
			dp.m(k,i,j).S.from = state_switch;
		} // k

	} // j
	} // i
}


mm_alignment::alignment_chain_ptr mm_alignment::viterbi_traceback( int k, int state, int i, int j,  dna_sequence_region_data &a, dna_sequence_region_data &b, uint16_t pos, bool chain_required, int sid ) {
	using namespace std;
	assert( pos <= 16384 );
  	/*
	cerr << "trace " << k << " ";
	switch( state ) {
		case state_M:
			cerr << "M " << dp.m(k,i,j).M.s.as_base();
			break;
		case state_Ia:
			cerr << "Ia " << dp.m(k,i,j).Ia.s.as_base();
			break;
		case state_Ib:
			cerr << "Ib " << dp.m(k,i,j).Ib.s.as_base();
			break;
	}
	cerr << endl;*/

	alignment_chain_ptr mychain(new alignment_chain());
	annotation_ptr empty_annotation( new annotation() );
	dna_alignment_sequence_ptr emptya( new dna_alignment_sequence()), emptyb( new dna_alignment_sequence() );
	mychain->ann = empty_annotation;
	mychain->align = pairwise_dna_alignment( emptya, emptyb, 0.0 );
	mychain->sid = sid;
	mychain->i_end = i;
	mychain->j_end = j;
	mychain->prev_i = 0;
	switch( state ) {
		case state_M:
			if( dp.m(k,i,j).M.s == 0.0 ) return mychain;
			break;
		case state_Ib:
			if( dp.m(k,i,j).Ib.s == 0.0 ) return mychain;
			break;
		case state_Ia:
			if( dp.m(k,i,j).Ia.s == 0.0 ) return mychain;
			break;
		case state_S:
			if( dp.m(k,i,j).S.s == 0.0 ) return mychain;
			break;
		case state_E:
			if( dp.m(k,i,j).E.s == 0.0 ) return mychain;
			break;
		case state_RIa:
			if( dp.r(i,j).Ia.s == 0.0 ) return mychain;
			break;
		case state_RIb:
			if( dp.r(i,j).Ib.s == 0.0 ) return mychain;
			break;
		case state_switch:
			if( dp.s(i,j).s == 0.0 ) return mychain;
	}

   // alignment sequences
   dna_alignment_sequence_ptr newa( new dna_alignment_sequence()), newb( new dna_alignment_sequence() );

	// annotation
	annotation_ptr myann( new annotation() ), final_ann( new annotation() );
	int rlast_col = -1, rlast_type = k;

   int next_i = i, next_j = j, next_k = k, next_state = state;
	int column = -1;
	bool chain = false, unannotated_columns = false;
	while( i > 0 || j > 0 || (!chain_required && (state != state_switch)) || (chain_required && !chain) ) {
		assert( i >= 0 );
		assert( j >= 0 );
      switch( state ) {
         case state_M:
            next_i = i - 1;
            next_j = j - 1;
				next_state = dp.m(k,i,j).M.from;
				if( (next_state >> 4) > 0 ) { chain = true; break; }
				assert( dp.m(k,i,j).M.from <= 16 );
				dp.m(k,i,j).M.from |= pos << 4;
            newa->data.push_front( a[(uint)i-1] );
            newb->data.push_front( b[(uint)j-1] );
				++column;
				unannotated_columns = true;
            break;
         case state_Ib:
            next_j = j - 1;
            next_state = dp.m(k,i,j).Ib.from;
				if( (next_state >> 4) > 0 ) { chain = true; break; }
				assert( dp.m(k,i,j).Ib.from <= 16 );
				dp.m(k,i,j).Ib.from |= pos << 4;
            newa->data.push_front( dna_alignment_alpha::GAP );
            newb->data.push_front( b[(uint)j-1] );
				++column;
				unannotated_columns = true;
            break;
         case state_Ia:
            next_i = i - 1;
            next_state = dp.m(k,i,j).Ia.from;
				if( (next_state >> 4) > 0 ) { chain = true; break; }
				assert( dp.m(k,i,j).Ia.from <= 16 );
				dp.m(k,i,j).Ia.from |= pos << 4;
            newa->data.push_front( a[(uint)i-1] );
            newb->data.push_front( dna_alignment_alpha::GAP );
				++column;
				unannotated_columns = true;
				break;
			case state_S:
				next_state = dp.m(k,i,j).S.from;
				break;
			case state_E:
				next_state = dp.m(k,i,j).E.from;
				break;
			case state_RIa:
				next_i = i - 1;
            next_state = dp.r(i,j).Ia.from;
				if( (next_state >> 4) > 0 ) { chain = true; break; }
				assert( dp.r(i,j).Ia.from <= 16 );
				dp.r(i,j).Ia.from |= pos << 4;
            newa->data.push_front( a[(uint)i-1] );
            newb->data.push_front( dna_alignment_alpha::GAP );
				++column;
				unannotated_columns = true;
				break;
			case state_RIb:
				next_j = j - 1;
            next_state = dp.r(i,j).Ib.from;
				if( (next_state >> 4) > 0 ) { chain = true; break; }
				assert( dp.r(i,j).Ib.from <= 16 );
				dp.r(i,j).Ib.from |= pos << 4;
            newa->data.push_front( dna_alignment_alpha::GAP );
            newb->data.push_front( b[(uint)j-1] );
				++column;
				unannotated_columns = true;
            break;
			case state_switch:
				if( dp.s(i,j).from == state_RIb ) {
					next_state = state_RIb;
					next_k = -1;
				} else {
					next_state = state_E;
					next_k = dp.s(i,j).from - state_E;
				}
				assert( unannotated_columns );

				// store annotated region
				if( rlast_type >= 0 ) {
					myann->push_back(annotation_interval(rlast_col+1,column,rlast_type+USER_MODEL));
				} else myann->push_back(annotation_interval(rlast_col+1,column,RANDOM_MODEL));
				unannotated_columns = false;

				//cerr << rlast_type << " " << next_k << " " << i << endl;
				rlast_col = column;
				rlast_type = next_k;
				break;

				// we encountered a previous chain
			default:
				cerr << state << "\t" << i << "\t" << j << endl;
				assert( 0 && "error in traceback" );
		}
		// preserve  values if we need to chain
		if( chain ) break;

      i = next_i;
      j = next_j;
		k = next_k;
		state = next_state;
   }
	if( chain_required && !chain ) cerr << "Did not chain: " << i << "\t" << j << "\t" << state << "\t" << k << endl;

	//cerr << " " << i << " " << j << " " << state << " " << k << endl;

	// store final region if needed
	if( unannotated_columns ) {
		if( rlast_type >= 0 ) {
			myann->push_back(annotation_interval(rlast_col+1,column,rlast_type+USER_MODEL));
		} else myann->push_back(annotation_interval(rlast_col+1,column,RANDOM_MODEL));
	}

	// now we need to invert the regions
	for( annotation_riter r = myann->rbegin(); r != myann->rend(); ++r ) {
		final_ann->push_back(annotation_interval(column - r->b, column - r->a, r->t));
	}

	mychain->align = pairwise_dna_alignment( newa, newb, 0.0 );
	mychain->ann = final_ann;
	mychain->i = i;
	mychain->j = j;
	mychain->prev_k = k;
	mychain->prev_i = (next_state >> 4);
	mychain->prev_state = state;

	return mychain;
}

void mm_alignment::do_statistics( std::pair<pairwise_dna_alignment,annotation_ptr> const &r, all_counts &cnt ) {
	assert( r.first.a->data.size() == r.first.b->data.size() );

	for( annotation_iter i = r.second->begin(); i != r.second->end(); ++i ) {

		// record region, here random is type RANDOM_MODEL and the models start at USER_MODEL
		cnt.region_uses[i->t] += 1;
		cnt.region_lengths[i->t] += i->b - i->a + 1;

		int gap_a = -1, gap_b = -1;
		for( int c = i->a; c <= i->b; ++c ) {
			assert( c < (int)r.first.a->data.size() );
			assert( c < (int)r.first.b->data.size() );


			dna_alignment_alpha::symbol c_a = r.first.a->data[c], c_b = r.first.b->data[c];

			// found a pair!
			if( c_a != dna_alignment_alpha::GAP && c_b != dna_alignment_alpha::GAP ) {
				int base_a = (int)dna_alpha::symbol(c_a);
				int base_b = (int)dna_alpha::symbol(c_b);
				cnt.pairs(i->t, base_a, base_b) += 1;
				cnt.bases_a(i->t,base_a) += 1;
				cnt.bases_b(i->t,base_b) += 1;

				// closing a gap in A
				if( gap_a >= 0 ) {
					cnt.gap_lengths_a[i->t] += c - gap_a + 1;
					gap_a = -1;
				}

				// closing a gap in B
				if( gap_b >= 0 ) {
					cnt.gap_lengths_b[i->t] += c - gap_b + 1;
					gap_b = -1;
				}

			// found a gap in A
			} else if( c_a == dna_alignment_alpha::GAP && c_b != dna_alignment_alpha::GAP ) {
				cnt.bases_b(i->t,(int)dna_alpha::symbol(c_b)) += 1;

				// gap is not open
				if( gap_a == -1 ) {
					gap_a = c;
					cnt.gaps_a[i->t] += 1;
				} // otherwise nothing to do

				// closing a gap in B
				if( gap_b >= 0 ) {
					cnt.gap_lengths_b[i->t] += c - gap_b + 1;
					gap_b = -1;
				}


			// found a gap in B
			} else if( c_a != dna_alignment_alpha::GAP && c_b == dna_alignment_alpha::GAP ) {
				cnt.bases_a(i->t,(int)dna_alpha::symbol(c_a)) += 1;

				// closing a gap in A
				if( gap_a >= 0 ) {
					cnt.gap_lengths_a[i->t] += c - gap_a + 1;
					gap_a = -1;
				}

				// gap is not open
				if( gap_b == -1 ) {
					gap_b = c;
				cnt.gaps_b[i->t] += 1;
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
			cnt.gap_lengths_a[i->t] += i->b - gap_a + 1;
			gap_a = -1;
		}

		// closing a gap in B
		if( gap_b >= 0 ) {
			cnt.gap_lengths_b[i->t] += i->b - gap_b + 1;
			gap_b = -1;
		}


	} // for i
}

void mm_alignment::perform_chain( mm_alignment::alignment_chain &chain, mm_alignment::interface_data &past, mm_alignment::interface_data &current ) {
	using namespace std;
	// must be starting region, nothing to chain
	if( chain.prev_i == 0 ) return;

	// if in interface, chain to previous
	if( chain.j == 0 && chain.i < interface*2 && past.alignments[chain.prev_i-1].get() != 0 ) chain.prev = past.alignments[chain.prev_i-1];
	else chain.prev = current.alignments[chain.prev_i-1];

}

std::pair<pairwise_dna_alignment,annotation_ptr> mm_alignment::collapse_chain( alignment_chain &c ) {
	int prev_i = c.prev_i, i = c.i, j = c.j, state = c.prev_state;
	alignment_chain_ptr prev = c.prev;

	using namespace std;
	/*std::cerr << "collapse_chain " << i << "\t" << j << std::endl;
	for( int id = 0; id < c.align.a->data.size(); ++id ) cerr << c.align.a->data[id];
	cerr << endl;
	for( int id = 0; id < c.align.b->data.size(); ++id ) cerr << c.align.b->data[id];
	cerr << endl;*/
	
	int current_segment = c.sid;
	while( prev_i > 0 ) {
		int col = prev->align.a->data.size()-1, ba = prev->i_end, bb = prev->j_end;
	
		/*cerr << ba << "\t" << bb << '\t' << i << '\t' << j << endl;
		for( int id = 0; id < prev->align.a->data.size(); ++id ) cerr << prev->align.a->data[id];
		cerr << endl;
		for( int id = 0; id < prev->align.b->data.size(); ++id ) cerr << prev->align.b->data[id];
		cerr << endl;*/

		// first find column to start at if we need to
		//if( !(j == 0 && i < interface*2) ) {
		if( current_segment == prev->sid ) {
			for( ; ba > i || bb > j; --col ) {
				assert( col >= 0 );
				if( prev->align.a->data[col] != dna_alignment_alpha::GAP ) --ba;
				if( prev->align.b->data[col] != dna_alignment_alpha::GAP ) --bb;
			}
			assert( ba == i );
			assert( bb == j );
		} // otherwise we just start at the end
//		cerr << "adjusted " << ba << '\t' << bb << endl;

		// add annotation
		for( annotation_riter k = prev->ann->rbegin(); k != prev->ann->rend(); ++k ) {
			// just add the annotation if not near the end column
			if( k->b <= col ) c.ann->push_front( *k );

			// check if we need to join the annotation intervals
			if( k->a <= col && k->b > col ) c.ann->push_front( annotation_interval(k->a,col,k->t) );

		}

		// add the required columns
		for( ; col >= 0; --col ) {
			c.align.a->data.push_front(prev->align.a->data[col]);
			c.align.b->data.push_front(prev->align.b->data[col]);
			if( prev->align.a->data[col] != dna_alignment_alpha::GAP ) --i;
			if( prev->align.b->data[col] != dna_alignment_alpha::GAP ) --j;
			//cerr << i << '\t' << prev->i << '\t' << j << '\t' << prev->j << endl;
		}

		prev_i = prev->prev_i;
		i = prev->i;
		j = prev->j;
		state = prev->prev_state;
		current_segment = prev->sid;
		prev = prev->prev;
	}
	
	// now fix annotation
	int c_last = -1;
	for( annotation_iter k = c.ann->begin(); k != c.ann->end(); ++k ) {
		assert( k->b >= 0 );
		assert( k->a >= 0 );
		assert( k->b >= k->a );
		k->b = c_last + 1 + k->b - k->a;
		k->a = c_last + 1;
		c_last = k->b;
	}

	// now merge regions that don't jump
	for( annotation_iter pk = c.ann->begin(), k = ++(c.ann->begin()); k != c.ann->end(); ) {
		if( pk->t == k-> t ) {
			pk->b = k->b;
			k = c.ann->erase(k);
		} else {
			++pk;
			++k;
		}
	}

	return make_pair(c.align,c.ann);
}

void mm_alignment::update_probabilities( mm_alignment::all_counts &counts ) {
	using namespace blitz;

	int total_regions = 0;
	// we have nmodels plus the random region
	for( int k = 0; k < (USER_MODEL+nmodels); ++k ) total_regions += counts.region_uses[k];

	// get combined base count for both sequences
	base_count all_bases(USER_MODEL+nmodels,dna_alpha::SIZE);
	all_bases = counts.bases_a + counts.bases_b;

	// compute probablity dist on regions
	for( int k = 0; k < nmodels; ++k ) {
		using namespace std;
		// region 0 is the random region
		p_sw[k] = (double)counts.region_uses[USER_MODEL+k]/(double)total_regions;

		// set model parameters if we ever use this region and have enough training data
		int total_bases = sum( all_bases(USER_MODEL+k,Range::all()) );

		if( counts.region_uses[USER_MODEL+k] > 0 && total_bases > 1000 ) {
			nw_model_parameters newp;

			double t = 1.0 - 1.0/(double)(counts.region_lengths[USER_MODEL+k]/counts.region_uses[USER_MODEL+k]);

			// compute p for region k
			int total_pairs = sum( counts.pairs(USER_MODEL+k,Range::all(),Range::all()) );
			for( int i = 0; i < dna_alpha::SIZE; ++i ) {
				for( int j = 0; j < dna_alpha::SIZE; ++j ) {
					newp.p(i,j) = (double) counts.pairs(USER_MODEL+k,i,j)/(double)total_pairs;
				} // j
			} // i


			for( int i = 0; i < dna_alpha::SIZE; ++i ) {
				newp.q(i) = (double) all_bases(USER_MODEL+k,i)/(double)total_bases;
			}

			// gap open
			int total_gaps = counts.gaps_a[USER_MODEL+k] + counts.gaps_b[USER_MODEL+k];
			newp.pr_open = (double)total_gaps/(double)total_bases;

			// mean gap length
			newp.mean_gap_length = (double)(counts.gap_lengths_a[USER_MODEL+k]+counts.gap_lengths_b[USER_MODEL+k])/(double)total_gaps;

			// set parameters
			mp[k].set(newp,t);
		}
	}

	// compute probability dist for random region
	p_rand = (double)counts.region_uses[RANDOM_MODEL]/(double)total_regions;

	// gaps in A are insertions in B
	if( counts.gaps_a[RANDOM_MODEL] > 5 ) {
		double mean_len_b = (double)counts.gap_lengths_a[RANDOM_MODEL]/(double)counts.gaps_a[RANDOM_MODEL];
		rp.p_b = 1.0 - 1.0/mean_len_b;
		rp.p_1mb = 1.0/mean_len_b;
	}

	// gaps in B are insertions in A
	if( counts.gaps_b[RANDOM_MODEL] > 5 ) {
		double mean_len_a = (double)counts.gap_lengths_b[RANDOM_MODEL]/(double)counts.gaps_b[RANDOM_MODEL];

		rp.p_a = 1.0 - 1.0/mean_len_a;
		rp.p_1ma = 1.0/mean_len_a;
	}

	// set model parameters if we ever use this region and have enough training data
	int total_bases = sum( all_bases(RANDOM_MODEL,Range::all()) );
	if( total_bases > 1000 ) {
		for( int i = 0; i < dna_alpha::SIZE; ++i ) {
			rp.q(i) = (double) all_bases((int)RANDOM_MODEL,i)/(double)total_bases;
		}
	}


}


void mm_alignment::set_interface_size( int i ) {
	assert( i > 0 );
	interface = i;
}

void mm_alignment::read_parameters( std::istream &in ) {
	using namespace std;
	vector<double> mdist;
	mdist.push_back( read<double>(in) );
	int rand_a = (int)read<double>(in);
	int rand_b = (int)read<double>(in);
	blitz::TinyVector<double,4> q;
	q = read<double>(in), read<double>(in), read<double>(in), read<double>(in);
	set_random_param(rand_a, rand_b, q);

	int nmodels = read<int>(in);
	for( int j = 0; j < nmodels; ++j ) {
		mdist.push_back( read<double>(in) );
		int length = (int)read<double>(in);
		nw_model_parameters param;
		in >> param;
		set_parameters( j, param, length );
	}
	set_region_dist( mdist );
}



void mm_alignment::start_training() {
	my_counts.reset( nmodels );
	my_pr = 0.0;
}

double mm_alignment::end_training() {
	update_probabilities( my_counts );
	return my_pr;
}

