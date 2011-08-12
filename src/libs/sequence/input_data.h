/*
 *  input_data.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-18.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef INPUT_DATA_H__
#define INPUT_DATA_H__
#include "dna_sequence.h"
#include "sequence_statistics.h"
#include "sequence_statistics_functions.h"

// information about a single input sequence
struct input_sequence_data {
	dna_sequence_ptr        seq;
	sequence_statistics_ptr stats;
};

// input data
struct input_data {
	input_sequence_data a;
	input_sequence_data b;
	input_data( dna_sequence_ptr ina, dna_sequence_ptr inb ) {
		a.seq = ina; b.seq = inb;
		a.stats = measure_statistics( *(a.seq) );
		b.stats = measure_statistics( *(b.seq) );
	}
};

#endif
