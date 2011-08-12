/*
 *  ex_sequence.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-09-16.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include <btl/sequence.h>
#include <btl/mdna_iupac_gap.h>
#include <boost/shared_ptr.hpp>
#include <vector>

typedef btl::sequence< btl::mdna_iupac_gap > ex_sequence;
typedef boost::shared_ptr< ex_sequence > ex_sequence_ptr;
typedef std::vector<ex_sequence_ptr> ex_sequence_ptr_vector;
typedef boost::shared_ptr< ex_sequence_ptr_vector > ex_sequence_ptr_vector_ptr;

