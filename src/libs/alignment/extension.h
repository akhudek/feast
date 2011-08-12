/*
 * extension.h
 *
 *  Created on: 29-Apr-2009
 *      Author: akhudek
 */

#ifndef EXTENSION_H_
#define EXTENSION_H_
#include <boost/shared_ptr.hpp>
#include "fragment.h"
#include "ereal.h"

struct extension {
	fragment_ptr_list 	segments;
	ereal 				score;
};

typedef boost::shared_ptr<extension> extension_ptr;

#endif /* EXTENSION_H_ */
