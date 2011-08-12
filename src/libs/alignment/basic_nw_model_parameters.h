#ifndef BASIC_NW_MODEL_PARAMETERS_H__
#define BASIC_NW_MODEL_PARAMETERS_H__
#include <blitz/tinyvec.h>

// This is mostly a data container, but it has attached methods.
struct basic_nw_model_parameters {	
	blitz::TinyVector<double,4> q;
   double pr_open;
   double mean_gap_length;
   
	// basic constructor
	basic_nw_model_parameters();
};

#endif
