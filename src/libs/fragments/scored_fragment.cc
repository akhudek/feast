/*
 *  scored_fragment.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-04-04.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "scored_fragment.h"

scored_fragment_ptr new_scored_fragment() {
	scored_fragment_ptr ptr( new scored_fragment() );
	return ptr;
}

