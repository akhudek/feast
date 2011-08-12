/*
 *  dynamic_fragment_score_worker.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#ifndef CUSTOM__FRAGMENT_SCORE_WORKER__
#define CUSTOM__FRAGMENT_SCORE_WORKER__
#include "nw_model_parameters.h"
#include "concurrent_queue.h"
#include "scored_fragment.h"
#include "input_data.h"

class custom_fragment_score_worker {
protected:
	std::vector<nw_model_parameters> parameter_sets;
	input_data &in;
	concurrent_queue<scored_fragment_ptr> &fragment_queue;
	
	local_segmentation<> lclseg;
	
	// segmentation extension 
	unsigned int EXT;
		
	// threshold at which to stop scoring a fragment
	ereal SHORT_CIRCUIT_THRESHOLD;
	
public:
	custom_fragment_score_worker( std::vector<nw_model_parameters> &ps, input_data &inseq, concurrent_queue<scored_fragment_ptr> &fq );
	
	void operator()();
	
};

#endif
