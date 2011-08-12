/*
 * diagonal_store.h
 *
 *  Created on: 18-Apr-2009
 *      Author: akhudek
 */

#ifndef DIAGONAL_STORE_H_
#define DIAGONAL_STORE_H_
#include <set>
#include "closed_interval.h"
#include "point.h"

class diagonal_store {
public:
	struct diagonal {
		int d;				// The diagonal.
		closed_interval r;	// Start and end on diagonal.
		diagonal( int nd, closed_interval nr ) : d(nd), r(nr) {}
	};

	struct diag_end_lt: public std::binary_function<diagonal, diagonal, bool> {
		diag_end_lt() {}

		inline bool operator()(diagonal const &I, diagonal const &J) const {
			if (I.d < J.d) return true;
			if (I.d == J.d && I.r.b < J.r.b) return true;
			return false;
		}
	};

private:
	typedef std::set<diagonal,diag_end_lt> store_t;
	store_t store;

public:
	diagonal_store();

	void add( int d, closed_interval nr );
	int compute_diagonal( point p ) const;
	bool overlaps_diagonal( point p ) const;
	void clear();
};



#endif /* DIAGONAL_STORE_H_ */
