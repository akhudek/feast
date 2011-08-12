/** Node and order functions for local chains.
 *
 */

#ifndef LOCAL_CHAINS_H_
#define LOCAL_CHAINS_H_

#include "ereal.h"

#include <btl/logspace.h>
#include <vector>
#include "point.h"


/** Node for a local chain. */
struct local_chain {
	/** Score of the chain ending here. */
	ereal score;

	/** Original point specifying a position in each of two sequences. */
	point p;

	/** Link to previous chain. */
	std::vector<local_chain>::iterator previous;

	/** Transformed point. */
	point tp;

	/** Skip pointer for faster DP in positive_local_chains algorithm. */
	std::vector<local_chain>::iterator xskip;

	/** For pruning chain tree. */
	bool seen;

	local_chain() : score(0.0), p(-1,-1), tp(-1,1 ), seen(false) {}
	local_chain( point &pp ) : score(0.0), p(pp), tp(pp), seen(false) {  tp.rotate(-3.14159265358979323846/4.0); }

};


struct local_chain_xylt: public std::binary_function<local_chain, local_chain, bool> {
	local_chain_xylt() {}

	inline bool operator()(local_chain const &I, local_chain const &J) {
		if (I.tp.x < J.tp.x) return true;
		if (I.tp.x == J.tp.x && I.tp.y < J.tp.y) return true;
		return false;
	}
};

struct local_chain_ylt: public std::binary_function<int, local_chain, bool> {
	local_chain_ylt() {}

	inline bool operator()(int const &v, local_chain const &J) {
		return J.tp.y < v;
	}
};

#endif /* LOCAL_CHAINS_H_ */
