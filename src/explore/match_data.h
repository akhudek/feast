/*
 *  match_data.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-09-16.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#include <vector>
#include <boost/shared_ptr.hpp>
#include <functional>
#include "ex_sequence.h"
#include "closed_interval.h"


typedef std::vector<double> match_data;
typedef boost::shared_ptr< match_data > match_data_ptr;

struct max_function : public std::binary_function<double, double, double> {
    double operator()(double x, double y) { return std::max(x, y); }
};

struct coverage : public std::binary_function<double, double, double> {
    double operator()(double x, double y) { return y + 1; }
};

template<typename ComposeFcn>
void fill_match_data(match_data &mdata, ex_sequence_ptr_vector_ptr aln, closed_interval r, int window_size, ComposeFcn fcn) {
	using namespace btl;
	using namespace std;
	using namespace blitz;
	ex_sequence::data_type &A = (*aln)[0]->data, &B = (*aln)[1]->data;
	unsigned int md_i = r.a;

	// initialize window
	int match_count = 0;
	for (unsigned int i = 0; i < (unsigned int)window_size; ++i) {
		if (A[i] != mdna_iupac_gap::GAP && B[i] != mdna_iupac_gap::GAP && A[i] == B[i]) match_count++;
	}

	assert( md_i < mdata.size() );
	mdata[md_i] = fcn((double) match_count / (double) window_size, mdata[md_i]);
	++md_i;

	for( unsigned int i = window_size, j = 0; i < A.size(); ++i, ++j) {
		if (A[j] != mdna_iupac_gap::GAP && B[j] != mdna_iupac_gap::GAP && A[j] == B[j]) 	match_count--;
		if (A[i] != mdna_iupac_gap::GAP && B[i] != mdna_iupac_gap::GAP && A[i] == B[i]) 	match_count++;

		assert( md_i < mdata.size() );
		mdata[md_i] = fcn((double) match_count / (double) window_size, mdata[md_i]);
		++md_i;
	}
}
