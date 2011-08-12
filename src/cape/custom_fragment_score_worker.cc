/*
 *  dynamic_fragment_score_worker.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "custom_fragment_score_worker.h"
#include "fragment_functions.h"
#include "fragment_algorithm.h"

custom_fragment_score_worker::custom_fragment_score_worker(
																			  std::vector<nw_model_parameters> &ps,
																			  input_data &indata,
																			  concurrent_queue<scored_fragment_ptr> &fq )
: parameter_sets(ps), in(indata), fragment_queue(fq), EXT(200)  {
	SHORT_CIRCUIT_THRESHOLD.set_base(200.0);
}

void custom_fragment_score_worker::operator()() {
	using namespace std;

	scored_fragment_ptr f;
	while( fragment_queue.get_data(f) ) {
		// skip if already scored
		if(f->score() > 0.0 ) continue;

		local_segmentation<>::result best_result;
		unsigned int best_parameter_set = 0;

		for( unsigned int j = 0; j < parameter_sets.size(); ++j ) {

			// estimate region composition
			scored_fragment newf = *f;

			lclseg.set_parameters( parameter_sets[j] );
			lclseg.set_q( parameter_sets[j].q );

         // compute forward and backwards extensions
         pair<int,int> midpoint = fragment_algorithm::midpoint( *f );
         local_segmentation<>::result fresult = lclseg.extend_forward( *in.a.seq, *in.b.seq, midpoint.first+1, midpoint.second+1, EXT );
         if( fresult.score > best_result.score ) {
				best_result = fresult;
				best_parameter_set = j;
			}
			//cerr << best_result.score << "\t" << j << endl;

			// obviously a real fragment so don't bother with anything more
			if( best_result.score > SHORT_CIRCUIT_THRESHOLD ) break;
      }

      // sort by score
		f->set_score_data( best_result );
      f->set_parameter_set( best_parameter_set );
	}
}
