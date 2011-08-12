#include "mm_trainer.h"
#include "ereal.h"
#include <blitz/tinyvec.h>
#include "io.h"
#include <iomanip>
#include <limits>

mm_trainer::mm_trainer( int i ) : nmodels(i), PHMM(kstates+i*SUBSIZE), dirty(true) {
	set_default_parameters( blitz::TinyVector<double,4>(0.25,0.25,0.25,0.25) );
}

void mm_trainer::set_default_parameters( blitz::TinyVector<double,4> q ) {
	using namespace std;
	using namespace boost;
	assert(nmodels < 30);

	typedef boost::tuple<double,double,double> double3;
	vector<double3> custom_pts;
	custom_pts.push_back( double3(0.4,0.0625,2.0) );
	custom_pts.push_back( double3(0.2,0.03125,1.333) );
	custom_pts.push_back( double3(0.6,0.0625,3.41424) );
	custom_pts.push_back( double3(0.8,0.03125,2.0) );

	// set random parameters
	double pr_rand = 1.0/10000.0;
	PHMM.set_transition(state_switch,state_RIa,pr_rand);
	set_random_param( 500, 500, q );


	assert( (unsigned int)nmodels <= custom_pts.size() );
	for( int k = 0; k < min(2,nmodels); ++k ) {
		nw_model_parameters newp;
		newp.q = q;
		newp.set_p_hky( custom_pts[k].get<0>(), 1.0 );
		newp.pr_open = custom_pts[k].get<1>();
		newp.mean_gap_length = custom_pts[k].get<2>();
		set_parameters( k, newp, 200 );

		int ms = kstates+k*SUBSIZE;
		// set switch
		PHMM.set_transition(state_switch,ms+state_S,(1.0-pr_rand)/nmodels);

		// set end to switch
		PHMM.set_transition(ms+state_E,state_switch,1.0);

		// set silent states for model 0
		PHMM.set_emmission(ms+state_S,APHMM::NOTHING,APHMM::NOTHING,1.0);
		PHMM.set_emmission(ms+state_E,APHMM::NOTHING,APHMM::NOTHING,1.0);
	}

	// model repeats
	for( int k = 2; k < min(5,nmodels); ++k ) {
		nw_model_parameters newp;
		switch(k) {
			case 2:
				newp.q = 0.4, 0.4, 0.1, 0.1;
				break;
			case 3:
				newp.q = 0.4, 0.1, 0.1, 0.4;
				break;
			case 4:
				newp.q = 0.1, 0.1, 0.4, 0.4;
				break;
		}
		newp.set_p_hky( 0.2, 1.0 );
		newp.pr_open = 0.03125;
		newp.mean_gap_length = 10;
		set_parameters( k, newp, 200 );

		int ms = kstates+k*SUBSIZE;
		// set switch
		PHMM.set_transition(state_switch,ms+state_S,(1.0-pr_rand)/nmodels);

		// set end to switch
		PHMM.set_transition(ms+state_E,state_switch,1.0);

		// set silent states for model 0
		PHMM.set_emmission(ms+state_S,APHMM::NOTHING,APHMM::NOTHING,1.0);
		PHMM.set_emmission(ms+state_E,APHMM::NOTHING,APHMM::NOTHING,1.0);
	}

	// set silent states for switch and random
	PHMM.set_emmission(state_switch,APHMM::NOTHING,APHMM::NOTHING,1.0);

	// set start and end states
	PHMM.set_start_state(state_switch);
	PHMM.set_end_state(state_switch);
	PHMM.optimize();
	dirty = false;
}

void mm_trainer::start_training() { PHMM.start_training(); }

ereal mm_trainer::end_training() { return PHMM.end_training(); }


int mm_trainer::num_models() const {
	return nmodels;
}

double mm_trainer::get_random_pr() const {
	return (double)PHMM.get_transition(state_switch,state_RIa);
}

double mm_trainer::get_random_length_a() const {
	return 1.0/(1.0 - (double)PHMM.get_transition(state_RIa,state_RIa));
}

double mm_trainer::get_random_length_b() const {
	return 1.0/(1.0 - (double)PHMM.get_transition(state_RIb,state_RIb));
}

double mm_trainer::get_region_pr( int k ) const {
	return (double)PHMM.get_transition(state_switch,kstates + k*SUBSIZE + state_S);
}

double mm_trainer::get_region_length( int k ) const {
	// since we didn't tie the "t" parts, we'll just average all three values
	int ms = kstates + k*SUBSIZE;
	double sum = (double)PHMM.get_transition(ms+state_M,ms+state_E);
	sum += (double)PHMM.get_transition(ms+state_Ia,ms+state_E);
	sum += (double)PHMM.get_transition(ms+state_Ib,ms+state_E);
	double p_1mt = sum/3.0;
	return 1.0/p_1mt;
}

void mm_trainer::set_random_param( int ua, int ub, blitz::TinyVector<double,4> q ) {
	dirty = true;
	double p_a = 1.0 - 1.0/(double)ua;
	double p_1ma = 1.0/(double)ua;
	double p_b = 1.0 - 1.0/(double)ub;
	double p_1mb = 1.0/(double)ub;
	PHMM.set_transition(state_RIa,state_RIa,p_a);
	PHMM.set_transition(state_RIa,state_RIb,p_1ma);
	PHMM.set_transition(state_RIb,state_RIb,p_b);
	PHMM.set_transition(state_RIb,state_switch,p_1mb);

	for( int i = 0; i < 4; ++i) {
		PHMM.set_emmission(state_RIa,i,APHMM::NOTHING,q(i));
		PHMM.set_emmission(state_RIb,APHMM::NOTHING,i,q(i));
	}
}

void mm_trainer::set_parameters( int k, nw_model_parameters const &nwp, int len ) {
	dirty = true;
	double pt = 1.0-1.0/(double)len;
	int ms = kstates + k*SUBSIZE;

	// set open
    double p_at = nwp.pr_open*pt;
	double p_1m2at = (1.0-2.0*nwp.pr_open)*pt;
	PHMM.set_transition(ms+state_M,ms+state_Ia,p_at);
	PHMM.set_transition(ms+state_M,ms+state_Ib,p_at);
	PHMM.set_transition(ms+state_M,ms+state_M,p_1m2at);

	// set start
	double p_a = nwp.pr_open;
	double p_1m2a = 1.0 - 2.0*nwp.pr_open;
	PHMM.set_transition(ms+state_S,ms+state_Ia,p_a);
	PHMM.set_transition(ms+state_S,ms+state_Ib,p_a);
	PHMM.set_transition(ms+state_S,ms+state_M,p_1m2a);

	// Since a gap is length 1 by the virtue of opening it, we adjust n down by 1.
	// The remaining gap extension probability corresponds to the mean of a
	// geometric distribution allowing zero extensions.
	double p_b = 1.0-1.0/nwp.mean_gap_length;
	double p_bt = p_b*pt;
	double p_1mbt = (1.0 - p_b)*pt;
	PHMM.set_transition(ms+state_Ia,ms+state_Ia,p_bt);
	PHMM.set_transition(ms+state_Ia,ms+state_M,p_1mbt);
	PHMM.set_transition(ms+state_Ib,ms+state_Ib,p_bt);
	PHMM.set_transition(ms+state_Ib,ms+state_M,p_1mbt);

	// set end transitions
	double p_1mt = 1.0 - pt;
	PHMM.set_transition( ms+state_M,ms+state_E,p_1mt);
	PHMM.set_transition( ms+state_Ia,ms+state_E,p_1mt);
	PHMM.set_transition( ms+state_Ib,ms+state_E,p_1mt);

	// set background
	for( int i = 0; i < 4; i++ ) {
      PHMM.set_emmission(ms+state_Ia,i,APHMM::NOTHING,nwp.q(i));
      PHMM.set_emmission(ms+state_Ib,APHMM::NOTHING,i,nwp.q(i));
   }

	// set p
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			PHMM.set_emmission(ms+state_M,i,j,nwp.p(i,j));
		}
   }
}

nw_model_parameters_ptr mm_trainer::get_parameters( int k ) const {
	nw_model_parameters_ptr nwp( new nw_model_parameters() );
	int ms = kstates + SUBSIZE*k;
	double t = 1.0 - 1.0/get_region_length(k);
	nwp->pr_open = (double) PHMM.get_transition(ms+state_M,ms+state_Ia)/t;
	nwp->pr_open += (double) PHMM.get_transition(ms+state_M,ms+state_Ib)/t;
	nwp->pr_open /= 2.0;

	double gap_a = 1.0/(1.0- (double)PHMM.get_transition(ms+state_Ia,ms+state_Ia)/t );
	double gap_b = 1.0/(1.0- (double)PHMM.get_transition(ms+state_Ib,ms+state_Ib)/t );
	nwp->mean_gap_length = (gap_a + gap_b)/2.0;

	for( int x = 0; x < APHMM::NOTHING; ++x ) {
		double qx = 0.0;
		for( int y = 0; y < APHMM::NOTHING; ++y ) {
			nwp->p(x,y) = (double)PHMM.get_emmission(ms+state_M,x,y);
			qx += nwp->p(x,y);
		}
		nwp->q(x) = qx;
	}

	return nwp;
}

blitz::TinyVector<double,4> mm_trainer::get_random_distribution() const {
	blitz::TinyVector<double,4> q;
	q = 0,0,0,0;
	for( int i = 0; i < 4; ++i ) {
		q(i) += (double)PHMM.get_emmission(state_RIa,i,APHMM::NOTHING);
		q(i) += (double)PHMM.get_emmission(state_RIb,APHMM::NOTHING,i);
		q(i) /= 2.0;
	}
	return q;
}

void mm_trainer::pretty_print_parameters( std::ostream &out ) const {
	using namespace std;

	out << "Random: " << get_random_pr() << endl;
	out << "Random length A: " << get_random_length_a() << " B: " << get_random_length_b() << endl;
	out  << "Random background: " << get_random_distribution() << endl;

	for( int k = 0; k < num_models(); ++k ) {
		out << endl << "Model " << k << ": " << get_region_pr(k) << endl;
		out << "Length: " << get_region_length(k) << endl;
		nw_model_parameters_ptr nwp = get_parameters(k);
		nwp->pretty_print(out);
	}
}

void mm_trainer::write_parameters( std::ostream &out ) const {
	using namespace std;

	out << "# Random (region 0) with pr " << get_random_pr();
	out << " and lengths A: " << get_random_length_a() << ", B: " << get_random_length_b() << endl;
	out  << "# Random background: " << get_random_distribution() << endl << endl;

	for( int j = 0; j < num_models(); ++j ) {
		out << "# Region " << (j+1) << " with pr " << get_region_pr(j);
		out << " and length " << get_region_length(j) << "." << endl;
		nw_model_parameters_ptr nwp = get_parameters(j);
		nwp->pretty_print( out, "# " );
		out << endl;
	}
	out << endl << "# Machine readable parameters follow. " << "---------------------------------------------" << endl;

	streamsize original_precision = out.precision();

	out << scientific << setprecision(numeric_limits<double>::digits10);

	out << get_random_pr() << " " << get_random_length_a() << " " << get_random_length_b() << " ";
	blitz::TinyVector<double,4> q = get_random_distribution();
	for( int i = 0; i < 4; ++i )  out << q(i) << " ";
	out << endl;

	out << num_models() << endl;
	for( int j = 0; j < num_models(); ++j ) {
		out <<  get_region_pr(j) << " " << get_region_length(j) << " " << *get_parameters(j) << endl;
	}

	out << resetiosflags(ios::scientific) << setprecision(original_precision);
}

void mm_trainer::read_parameters( std::istream &in ) {
	using namespace std;
	vector<double> mdist;
	mdist.push_back( read<double>(in) );
	int rand_a = (int)read<double>(in);
	int rand_b = (int)read<double>(in);
	blitz::TinyVector<double,4> q;
	for( int i = 0; i < 4; ++i ) q(i) = read<double>(in);
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
	PHMM.optimize();
	dirty = false;
}

void mm_trainer::disable_random_model() {
	// we redistribution the probability assigned to the random model to the remaining
	// models evenly
	double prslice = (double)PHMM.get_transition(state_switch,state_RIa)/(double)nmodels;
	for( int k = 0; k < nmodels; ++k ) {
		double pr = PHMM.get_transition(state_switch,kstates + k*SUBSIZE);
		PHMM.set_transition(state_switch,kstates + k*SUBSIZE, pr + prslice);
	}
	PHMM.set_transition(state_switch,state_RIa, 0.0 );

	// set all these to zero to allow our optimizer to ignore them
	PHMM.set_transition(state_RIa,state_RIa, 0.0 );
	PHMM.set_transition(state_RIa,state_RIb, 0.0 );
	PHMM.set_transition(state_RIb,state_RIb, 0.0 );
	PHMM.set_transition(state_RIb,state_switch, 0.0 );

	PHMM.optimize();
	dirty = false;
}
