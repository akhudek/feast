/*
 * scored_fragment_selectors.h
 *
 *  Created on: 12-Dec-2008
 *      Author: akhudek
 */

#ifndef SCORED_FRAGMENT_SELECTORS_H_
#define SCORED_FRAGMENT_SELECTORS_H_

#include "scored_fragment.h"

struct by_increasing_Bstart : public std::binary_function<scored_fragment_ptr,scored_fragment_ptr,bool> {
public:
	inline bool operator() ( scored_fragment_ptr a, scored_fragment_ptr b ) {
		return (*a)[1].a < (*b)[1].a;
	}
};


#endif /* SCORED_FRAGMENT_SELECTORS_H_ */
