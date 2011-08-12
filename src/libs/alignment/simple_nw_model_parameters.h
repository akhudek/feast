/*
 *  simple_nw_model_parameters.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#ifndef SIMPLE_NW_MODEL_PARAMETERS_H__
#define SIMPLE_NW_MODEL_PARAMETERS_H__
#include <blitz/tinyvec.h>
#include "basic_nw_model_parameters.h"

// This is mostly a data container, but it has attached methods.
struct simple_nw_model_parameters : public basic_nw_model_parameters {
	double d;
   
	// basic constructor
	simple_nw_model_parameters();
	
	simple_nw_model_parameters( double d, double popen, double mgl );

};


#endif
