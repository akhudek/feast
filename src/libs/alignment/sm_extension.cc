#include "sm_extension.h"
#include <blitz/tinyvec.h>
#include "io.h"

sm_extension::sm_extension() : PHMM(NUM_STATES), dirty(true), stop_at_repeat(false) {
	set_default_parameters( blitz::TinyVector<double,4>(0.25,0.25,0.25,0.25) );
}

// Legacy to create a uniform api with mm_alignment.
void sm_extension::disable_random_model() {}

int sm_extension::num_models() const { return 1; }

void sm_extension::reset() { store.clear(); dstore.clear(); }

void sm_extension::set_default_parameters( blitz::TinyVector<double,4> nq ) {
	using namespace std;
	using namespace boost;

	nw_model_parameters newp;
	newp.q = nq;
	newp.set_p_hky( 0.4, 1.0 );
	newp.pr_open = 0.0625;
	newp.mean_gap_length = 2.0;
	set_parameters( 0, newp, 200 );

	// set start and end states
	PHMM.set_start_state(state_M);
	PHMM.set_end_state(state_E);

	// set silent states for end state
	PHMM.set_emmission(state_E,APHMM::NOTHING,APHMM::NOTHING,1.0);

	set_random_distribution(nq);
	PHMM.optimize();
	dirty = false;
}

void sm_extension::set_stop_at_repeat( bool t ) { stop_at_repeat = t; }

void sm_extension::set_random_distribution( blitz::TinyVector<double,4> nq ) {
	dirty = true;
	for( int i = 0; i < 4; ++i) q(i) = nq(i);
}

double sm_extension::get_region_length( int k ) const {
	assert( k == 0 );
	// since we didn't tie the "t" parts, we'll just average all three values
	double sum = (double)PHMM.get_transition(state_M,state_E);
	sum += (double)PHMM.get_transition(state_Ia,state_E);
	sum += (double)PHMM.get_transition(state_Ib,state_E);
	double p_1mt = sum/3.0;
	return 1.0/p_1mt;
}

void sm_extension::set_parameters( int k, nw_model_parameters const &nwp, int len ) {
	assert( k == 0 );
	dirty = true;
	double pt = 1.0-1.0/(double)len;

	// set open
	double p_at = nwp.pr_open*pt;
	double p_1m2at = (1.0-2.0*nwp.pr_open)*pt;
	PHMM.set_transition(state_M,state_Ia,p_at);
	PHMM.set_transition(state_M,state_Ib,p_at);
	PHMM.set_transition(state_M,state_M,p_1m2at);

	// Since a gap is length 1 by the virtue of opening it, we adjust n down by 1.
	// The remaining gap extension probability corresponds to the mean of a
	// geometric distribution allowing zero extensions.
	double p_b = 1.0-1.0/nwp.mean_gap_length;
	double p_bt = p_b*pt;
	double p_1mbt = (1.0 - p_b)*pt;
	PHMM.set_transition(state_Ia,state_Ia,p_bt);
	PHMM.set_transition(state_Ia,state_M,p_1mbt);
	PHMM.set_transition(state_Ib,state_Ib,p_bt);
	PHMM.set_transition(state_Ib,state_M,p_1mbt);

	// set end transitions
	double p_1mt = 1.0 - pt;
	PHMM.set_transition(state_M,state_E,p_1mt);
	PHMM.set_transition(state_Ia,state_E,p_1mt);
	PHMM.set_transition(state_Ib,state_E,p_1mt);

	// set background
	for( int i = 0; i < 4; i++ ) {
      PHMM.set_emmission(state_Ia,i,APHMM::NOTHING,nwp.q(i));
      PHMM.set_emmission(state_Ib,APHMM::NOTHING,i,nwp.q(i));
   }

	// set p
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			PHMM.set_emmission(state_M,i,j,nwp.p(i,j));
		}
   }
}

nw_model_parameters_ptr sm_extension::get_parameters( int k ) const {
	assert( k == 0 );
	nw_model_parameters_ptr nwp( new nw_model_parameters() );
	double t = 1.0 - 1.0/get_region_length(0);
	nwp->pr_open = (double) PHMM.get_transition(state_M,state_Ia)/t;

	nwp->mean_gap_length = 1.0/(1.0- (double)PHMM.get_transition(state_Ia,state_Ia)/t );

	for( int x = 0; x < APHMM::NOTHING; ++x ) {
		double qx = 0.0;
		for( int y = 0; y < APHMM::NOTHING; ++y ) {
			nwp->p(x,y) = (double)PHMM.get_emmission(state_M,x,y);
			qx += nwp->p(x,y);
		}
		nwp->q(x) = qx;
	}

	return nwp;
}

blitz::TinyVector<double,4> sm_extension::get_random_distribution() const {
	return q;
}

void sm_extension::pretty_print_parameters( std::ostream &out ) const {
	using namespace std;

	out  << "Background: " << get_random_distribution() << endl;
	nw_model_parameters_ptr nwp = get_parameters(0);
	nwp->pretty_print(out);
	double pr_match = 0;
	for( int i = 0; i < 4; ++i ) pr_match += nwp->p(i,i);
	out << "Pr(match): " << pr_match << endl;
}

void sm_extension::write_parameters( std::ostream &out ) const {
	using namespace std;

	out  << "# Background: " << get_random_distribution() << endl << endl;

	nw_model_parameters_ptr nwp = get_parameters(0);
	nwp->pretty_print( out, "# " );
	out << endl;
	out << endl << "# Machine readable parameters follow. " << "---------------------------------------------" << endl;

	blitz::TinyVector<double,4> q = get_random_distribution();
	for( int i = 0; i < 4; ++i )  out << q(i) << " ";
	out << endl;
	out << *get_parameters(0) << endl;

}

void sm_extension::read_parameters( std::istream &in ) {
	using namespace std;
	vector<double> mdist;
	mdist.push_back( read<double>(in) );
	// Skip unused random model parameters.
	read<double>(in);	read<double>(in);
	for( int i = 0; i < 4; ++i ) q(i) = read<double>(in);

	int nmodels = read<int>(in);
	assert(nmodels == 1);

	mdist.push_back( read<double>(in) );
	int length =  read<double>(in);
	nw_model_parameters param;
	in >> param;
	set_parameters( 0, param, length );
	PHMM.optimize();
	dirty = false;
}
