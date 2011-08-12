/*
 *  basic_nw_model_parameters.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "basic_nw_model_parameters.h"

basic_nw_model_parameters::basic_nw_model_parameters() : pr_open(0.0), mean_gap_length(0.0) {
	q = 0.25, 0.25, 0.25, 0.25;
}
