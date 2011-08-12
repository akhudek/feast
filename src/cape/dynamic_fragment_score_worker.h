/*
 *  dynamic_fragment_score_worker.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#ifndef DYNAMIC_FRAGMENT_SCORE_WORKER__
#define DYNAMIC_FRAGMENT_SCORE_WORKER__
#include "simple_nw_model_parameters.h"
#include "concurrent_queue.h"
#include "scored_fragment.h"
#include "input_data.h"

class dynamic_fragment_score_worker {
protected:
	std::vector<simple_nw_model_parameters> parameter_sets;
	input_data &in;
	concurrent_queue<scored_fragment_ptr> &fragment_queue;

	local_segmentation<> lclseg;
	
	// segmentation extension 
	unsigned int EXT;
	
	// the minimum sized window we use to measure composition
	unsigned int min_composition_window;
	
	// threshold at which to stop scoring a fragment
	ereal SHORT_CIRCUIT_THRESHOLD;

public:
	dynamic_fragment_score_worker( std::vector<simple_nw_model_parameters> &ps, input_data &inseq, concurrent_queue<scored_fragment_ptr> &fq, bool use_mask );

	void operator()();
	
};

#endif
