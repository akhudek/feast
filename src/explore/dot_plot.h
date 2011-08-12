/*
 * dot_plot.h
 *
 *  Created on: 11-Feb-2009
 *      Author: akhudek
 */

#ifndef DOT_PLOT_H_
#define DOT_PLOT_H_

#include <ostream>
#include "ex_sequence.h"

void dot_plot_coords( std::ostream &out, ex_sequence_ptr_vector_ptr aln, int start_a, int start_b, int id );

#endif /* DOT_PLOT_H_ */
