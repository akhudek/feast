#include "nw_model_parameters.h"
#include "io.h"
#include "dna_sequence.h"

nw_model_parameters::nw_model_parameters() {
	p.resize(dna_alpha::alphabet.size(),dna_alpha::alphabet.size());
	q = 0.25, 0.25, 0.25, 0.25;
	mean_gap_length = 2.0;
	pr_open = 0.0625;
	set_p_hky( 0.4, 1.0 );
}

nw_model_parameters::nw_model_parameters( nw_model_parameters const &other ) {
	p.resize(other.p.extent());
	p = other.p;
	q = other.q;
	mean_gap_length = other.mean_gap_length;
	pr_open = other.pr_open;
}

nw_model_parameters::nw_model_parameters( simple_nw_model_parameters const &other, double ratio ) {
	p.resize(dna_alpha::alphabet.size(),dna_alpha::alphabet.size());
	q = other.q;
	mean_gap_length = other.mean_gap_length;
	pr_open = other.pr_open;
	set_p_hky( other.d, ratio );
}

nw_model_parameters &nw_model_parameters::operator=( nw_model_parameters const &other ) {
	p = other.p;
	q = other.q;
	mean_gap_length = other.mean_gap_length;
	pr_open = other.pr_open;
	return *this;
}

void nw_model_parameters::set_p_hky( double d, double a ) {
   using namespace std;

   // the default transition matrix needs to be based on q
   // NOTE: I did not type the following out. I am not crazy. This was copy and pasted from maple to ensure correctness.
   p((int)dna_alpha::A,(int)dna_alpha::A) = (-1.0+2.0*q((int)dna_alpha::C)+q((int)dna_alpha::G)+2.0*q((int)dna_alpha::T)
															-pow(q((int)dna_alpha::T),2.0)-pow(q((int)dna_alpha::C),2.0)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
															*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																			  +q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																			  *pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															-q((int)dna_alpha::T)*exp(0.5*a*d/(-q((int)dna_alpha::G)
																										  +q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a
																										  *q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																										  -q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))+2.0*exp(0.5*a
																																																														  *d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																																																+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																																																																*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)
																																																																*a))*q((int)dna_alpha::T)*q((int)dna_alpha::C)-2.0*q((int)dna_alpha::T)*q((int)dna_alpha::C)+exp(0.5*a*d
																																																																																																 /(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)
																																																																																																	-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																																																																																																	-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))*pow(q((int)dna_alpha::C),2.0)
															-q((int)dna_alpha::G)*q((int)dna_alpha::T)-q((int)dna_alpha::C)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																																								  *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a
																																								  *q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																																								  -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-q((int)dna_alpha::C)*q((int)dna_alpha::G)+exp(0.5*a
																																																																											 *d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																																																																												  *q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																																																																												  -q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))*pow(q((int)dna_alpha::T),2.0)
															+q((int)dna_alpha::C)*q((int)dna_alpha::G)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																																	  *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																																	  +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)
																																	  *a+pow(q((int)dna_alpha::T),2.0)*a))-exp(0.5*(1-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a*q((int)dna_alpha::C)
																																																	+q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																																										 +q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																																																										 +a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															*q((int)dna_alpha::G))/(-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::C));

   p((int)dna_alpha::A,(int)dna_alpha::C) = -(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																					 +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																					 +a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																					 +pow(q((int)dna_alpha::T),2.0)*a)))*q((int)dna_alpha::C);

   p((int)dna_alpha::A,(int)dna_alpha::G) = -q((int)dna_alpha::G)*(1+q((int)dna_alpha::T)*exp(0.5*a*d/(-q((int)dna_alpha::G)
																																		 +q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																																		 +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																																		 -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-q((int)dna_alpha::T)-q((int)dna_alpha::C)+q((int)dna_alpha::C)
																						 *exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																											*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																											*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)
																											*a))-exp(0.5*(1-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a*q((int)dna_alpha::C)+q((int)dna_alpha::T)*a)
																														*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																															 *q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																															 -q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)))/(-1.0+q((int)dna_alpha::T)
																																																																	 +q((int)dna_alpha::C));

   p((int)dna_alpha::A,(int)dna_alpha::T) = -(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																					 +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)
																					 *q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)
																					 *a+pow(q((int)dna_alpha::T),2.0)*a)))*q((int)dna_alpha::T);

   p((int)dna_alpha::C,(int)dna_alpha::A) = (-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::G)+q((int)dna_alpha::C))*(-1.0+exp(0.5*a
																																									 *d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																																										  *q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																																										  -q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)));

	p((int)dna_alpha::C,(int)dna_alpha::C) = -( -pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::C)
															 *exp(
																0.5*a*d/(
																	-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																	+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																	*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)
																	*a)
																)
															 +exp(
																0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)
																*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)
																*a+pow(q((int)dna_alpha::T),2.0)*a))
															 *q((int)dna_alpha::T)*q((int)dna_alpha::C)
															 +exp(
																	0.5*a*d/(
																	-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																	*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																	-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a
																	)
																)
															*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*exp(
																	-0.5*(
																		-q((int)dna_alpha::C)-a+q((int)dna_alpha::T)*a+a
																		*q((int)dna_alpha::C)-q((int)dna_alpha::T)
																		)
																	*d/(
																		-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																		+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)
																		*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																		+pow(q((int)dna_alpha::T),2.0)*a
																		)
																	)
																)
															 /(q((int)dna_alpha::C)+q((int)dna_alpha::T));

	p((int)dna_alpha::C,(int)dna_alpha::G) = -q((int)dna_alpha::G)*(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																												 *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																												 +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																												 -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)));

   p((int)dna_alpha::C,(int)dna_alpha::T) = -q((int)dna_alpha::T)*(-q((int)dna_alpha::C)-exp(0.5*a*d/(-q((int)dna_alpha::G)
																																		+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)
																																		-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)
																																		*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))+q((int)dna_alpha::T)
																						 *exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																											*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																											+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
																						 +q((int)dna_alpha::C)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																																		*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a
																																		*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																																		-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-q((int)dna_alpha::T)+exp(-0.5*(-q((int)dna_alpha::C)-a
																																																															  +q((int)dna_alpha::T)*a+a*q((int)dna_alpha::C)-q((int)dna_alpha::T))*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																																																																																							  *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																																																																																							  +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																																																																																							  -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)))/(q((int)dna_alpha::C)+q((int)dna_alpha::T));

   p((int)dna_alpha::G,(int)dna_alpha::A) = (-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::G)+q((int)dna_alpha::C))
	*(1+q((int)dna_alpha::T)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
													  +q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
													  +a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
	  -q((int)dna_alpha::T)-q((int)dna_alpha::C)+q((int)dna_alpha::C)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																											 *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a
																											 *q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																											 -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-exp(0.5*(1-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a
																																																	  *q((int)dna_alpha::C)+q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																																																																		+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																																																																		+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																																																																		+pow(q((int)dna_alpha::T),2.0)*a)))/(-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::C));

   p((int)dna_alpha::G,(int)dna_alpha::C) = -(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																					 +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)
																					 *q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																					 +pow(q((int)dna_alpha::T),2.0)*a)))*q((int)dna_alpha::C);

   p((int)dna_alpha::G,(int)dna_alpha::G) = -(q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)*exp(0.5*a*d
																																					  /(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																																						 *q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)
																																						 -q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															 -q((int)dna_alpha::G)*q((int)dna_alpha::T)+exp(0.5*(1.0-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a*q((int)dna_alpha::C)
																																  +q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																										+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																																										*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															 -q((int)dna_alpha::C)*q((int)dna_alpha::G)+q((int)dna_alpha::C)*q((int)dna_alpha::G)*exp(0.5*a*d/(-q((int)dna_alpha::G)
																																																+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																																																+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																																																-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-exp(0.5*(1.0-q((int)dna_alpha::T)-q((int)dna_alpha::C)
																																																																						 +a*q((int)dna_alpha::C)+q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																																																																																							 +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																																																																																							 +a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																																																																																							 +pow(q((int)dna_alpha::T),2.0)*a))*q((int)dna_alpha::C)-exp(0.5*(1.0-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a*q((int)dna_alpha::C)
																																																																																																													+q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																																																																																																						 +q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																																																																																																																						 *pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															 *q((int)dna_alpha::T)-exp(0.5*(1-q((int)dna_alpha::T)-q((int)dna_alpha::C)+a*q((int)dna_alpha::C)
																									  +q((int)dna_alpha::T)*a)*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																			+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)
																																			+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))
															 *q((int)dna_alpha::G))/(-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::C));

   p((int)dna_alpha::G,(int)dna_alpha::T) = -(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																					 +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a
																					 *q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																					 -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)))*q((int)dna_alpha::T);

   p((int)dna_alpha::T,(int)dna_alpha::A) = (-1.0+q((int)dna_alpha::T)+q((int)dna_alpha::G)+q((int)dna_alpha::C))*(-1.0+exp(0.5*a
																																									 *d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																																										  +q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																																										  *pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)));

   p((int)dna_alpha::T,(int)dna_alpha::C) = -q((int)dna_alpha::C)*(-q((int)dna_alpha::C)-exp(0.5*a*d/(-q((int)dna_alpha::G)
																																		+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a
																																		*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)
																																		*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))+q((int)dna_alpha::T)
																						 *exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																											+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																											*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)
																											*a))+q((int)dna_alpha::C)*exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																																								*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a
																																								*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)
																																								*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-q((int)dna_alpha::T)+exp(-0.5
																																																																										  *(-q((int)dna_alpha::C)-a+q((int)dna_alpha::T)*a+a*q((int)dna_alpha::C)-q((int)dna_alpha::T))
																																																																										  *d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)
																																																																												*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																																																																												*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)
																																																																												*a)))/(q((int)dna_alpha::C)+q((int)dna_alpha::T));

   p((int)dna_alpha::T,(int)dna_alpha::G) = -q((int)dna_alpha::G)*(-1.0+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)
																												 *q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																												 +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)
																												 -q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a)));

   p((int)dna_alpha::T,(int)dna_alpha::T) = -(-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::C)*exp(-0.5
																																						*(-q((int)dna_alpha::C)-a+q((int)dna_alpha::T)*a+a*q((int)dna_alpha::C)-q((int)dna_alpha::T))*d/(-q((int)dna_alpha::G)
																																																																						 +q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)
																																																																						 -a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)
																																																																						 *q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))-pow(q((int)dna_alpha::T),2.0)-q((int)dna_alpha::T)
															 *exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																				+pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)
																				*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)
																				*a+pow(q((int)dna_alpha::T),2.0)*a))+exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)
																																				  +pow(q((int)dna_alpha::G),2.0)+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)
																																				  +2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)
																																				  *q((int)dna_alpha::C)-q((int)dna_alpha::T)*a+pow(q((int)dna_alpha::T),2.0)*a))*pow(q((int)dna_alpha::T),2.0)
															 +exp(0.5*a*d/(-q((int)dna_alpha::G)+q((int)dna_alpha::G)*q((int)dna_alpha::T)+pow(q((int)dna_alpha::G),2.0)
																				+q((int)dna_alpha::C)*q((int)dna_alpha::G)-a*q((int)dna_alpha::C)+2.0*a*q((int)dna_alpha::C)*q((int)dna_alpha::T)+a
																				*pow(q((int)dna_alpha::C),2.0)-q((int)dna_alpha::T)*q((int)dna_alpha::C)-q((int)dna_alpha::T)*a
																				+pow(q((int)dna_alpha::T),2.0)*a))*q((int)dna_alpha::T)*q((int)dna_alpha::C))/(q((int)dna_alpha::C)+q((int)dna_alpha::T));

   p((int)dna_alpha::A,(int)dna_alpha::A) *= q((int)dna_alpha::A);
   p((int)dna_alpha::A,(int)dna_alpha::T) *= q((int)dna_alpha::A);
   p((int)dna_alpha::A,(int)dna_alpha::C) *= q((int)dna_alpha::A);
   p((int)dna_alpha::A,(int)dna_alpha::G) *= q((int)dna_alpha::A);

   p((int)dna_alpha::T,(int)dna_alpha::A) *= q((int)dna_alpha::T);
   p((int)dna_alpha::T,(int)dna_alpha::T) *= q((int)dna_alpha::T);
   p((int)dna_alpha::T,(int)dna_alpha::C) *= q((int)dna_alpha::T);
   p((int)dna_alpha::T,(int)dna_alpha::G) *= q((int)dna_alpha::T);

   p((int)dna_alpha::C,(int)dna_alpha::A) *= q((int)dna_alpha::C);
   p((int)dna_alpha::C,(int)dna_alpha::T) *= q((int)dna_alpha::C);
   p((int)dna_alpha::C,(int)dna_alpha::C) *= q((int)dna_alpha::C);
   p((int)dna_alpha::C,(int)dna_alpha::G) *= q((int)dna_alpha::C);

   p((int)dna_alpha::G,(int)dna_alpha::A) *= q((int)dna_alpha::G);
   p((int)dna_alpha::G,(int)dna_alpha::T) *= q((int)dna_alpha::G);
   p((int)dna_alpha::G,(int)dna_alpha::C) *= q((int)dna_alpha::G);
   p((int)dna_alpha::G,(int)dna_alpha::G) *= q((int)dna_alpha::G);
}

std::istream &operator>>( std::istream &in, nw_model_parameters &nw ) {
	// read p
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			nw.p(i,j) = read<double>(in);
		}
   }

	// read q
	for( int i = 0; i < 4; i++ ) {
      nw.q(i) = read<double>(in);
   }

	nw.pr_open = read<double>( in );
   nw.mean_gap_length = read<double>( in );
	return in;
}

void nw_model_parameters::pretty_print( std::ostream &out,  char const *prefix ) const {
	using namespace std;
	streamsize old_prec = out.precision();
	out << setprecision(4);
	out << prefix << "Substitution matrix (A T C G):" << endl;
	for( int i = 0; i < 4; i++ ) {
		out << prefix;
		for( int j = 0; j < 4; j++ ) {
			out << p(i,j) << " ";
		}
		out << endl;
   }
	out << prefix << "Background (A T C G): ";
	for( int i = 0; i < 4; ++i ) out << q(i) <<  " ";
	out << endl;
	out << prefix << "pr_open: " << pr_open << " cost: " << log(pr_open)/log(2.0) << " bits" << endl;
	out << prefix << "mean gap length: " << mean_gap_length;
	double ext = 1.0-1.0/mean_gap_length;
	out << " extend cost: " << log(ext)/log(2.0) << " bits" << endl;

	double pr_match = 0;
	for( int i = 0; i < 4; ++i ) pr_match += p(i,i);
	out << prefix << "Expected identities: " << (pr_match*100) << "%" << endl;
	out << setprecision(old_prec);
}

std::ostream &operator<<( std::ostream &out, nw_model_parameters const &nw ) {
	using namespace std;
	// read p
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			out << nw.p(i,j)  << " ";
		}
   }

	// read q
	for( int i = 0; i < 4; i++ ) {
      out << nw.q(i) << " ";
   }

	out << nw.pr_open << " ";
	out << nw.mean_gap_length << endl;
	return out;
}

