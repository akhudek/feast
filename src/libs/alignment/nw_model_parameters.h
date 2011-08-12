#ifndef NW_MODEL_PARAMETERS_H__
#define NW_MODEL_PARAMETERS_H__
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "basic_nw_model_parameters.h"
#include "simple_nw_model_parameters.h"

// This is mostly a data container, but it has attached methods.
struct nw_model_parameters : basic_nw_model_parameters {
	friend std::istream &operator>>( std::istream &in, nw_model_parameters &nw ); 
	friend std::ostream &operator<<( std::ostream &out, nw_model_parameters const &nw ); 
	
   blitz::Array<double,2> p;
	
	// basic constructor
	nw_model_parameters();
	
	// copy constructor
	nw_model_parameters( nw_model_parameters const &other );
	
	// copy constructor from simple
	nw_model_parameters( simple_nw_model_parameters const &other, double ratio = 1.0 );
	
	// copy operation
	nw_model_parameters &operator=( nw_model_parameters const &other );
	
	// Sets p using HKY with d subs/site and 'a' transversion/transition ratio.
	void set_p_hky( double d, double a );
	
	void pretty_print( std::ostream &out, char const *prefix = "" ) const;

};

typedef boost::shared_ptr<nw_model_parameters> nw_model_parameters_ptr;


#endif
