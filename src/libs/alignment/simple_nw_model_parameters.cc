/*
 *  simple_nw_model_parameters.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#include "simple_nw_model_parameters.h"


simple_nw_model_parameters::simple_nw_model_parameters() : basic_nw_model_parameters(), d(0.0) {}

simple_nw_model_parameters::simple_nw_model_parameters( double nd, double npo, double mgl ) : basic_nw_model_parameters(), d(nd) {
	pr_open = npo;
	mean_gap_length = mgl;
}

