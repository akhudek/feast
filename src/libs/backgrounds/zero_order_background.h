#ifndef ZERO_ORDER_BACKGROUND_H__
#define ZERO_ORDER_BACKGROUND_H__
#include "ereal.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <blitz/tinyvec.h>
#include <btl/logspace.h>
#include "closed_interval.h"
#include "dna_sequence.h"

class zero_order_background {
protected:
	// background probabilities over A,T,C,G
	blitz::TinyVector<ereal, 4> q;

	// structure storing preprocess data
	std::vector< ereal > data;

public:
	template<typename SEQ>
	zero_order_background(blitz::TinyVector<double, 4> newq, SEQ &a );

	template<typename SEQ>
	zero_order_background(blitz::TinyVector<double, 4> newq, SEQ &a, closed_interval r );


	// compute model probability for a subinterval given preprocessing data
	inline ereal pr( closed_interval r );
};

template<typename SEQ>
zero_order_background::zero_order_background( blitz::TinyVector<double, 4> newq, SEQ &a, closed_interval r ) {
	// set q
	for( int i = 0; i < 4; i++ ) q(i) = newq(i);

	// create data structure, reserve memory, and save the range
	data.reserve(a.size());

	// start the computation
	ereal last = 1.0;
	for( typename SEQ::const_iterator i = a.begin()+r.a; i != a.begin()+r.b+1; ++i ) {
		last *= q((int) *i);
		data.push_back(last);
	}
}


template<typename SEQ>
zero_order_background::zero_order_background( blitz::TinyVector<double, 4> newq, SEQ &a ) {
	// set q
	for( int i = 0; i < 4; i++ ) q(i) = newq(i);

	// create data structure, reserve memory, and save the range
	data.reserve(a.size());

	// start the computation
	ereal last = 1.0;
	for( typename SEQ::const_iterator i = a.begin(); i != a.end(); ++i ) {
		last *= q((int) *i);
		data.push_back(last);
	}
}



inline ereal zero_order_background::pr( closed_interval r) {
	if (r.b < r.a) return 1.0;
	return (r.a > 0) ? data[r.b] / data[r.a - 1] : data[r.b];
}

#endif
