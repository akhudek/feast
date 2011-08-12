#include "diagonal_store.h"

diagonal_store::diagonal_store() {}

void diagonal_store::add( int d, closed_interval nr ) {
	store.insert(diagonal(d,nr));
}

bool diagonal_store::overlaps_diagonal( point p ) const {
	using namespace std;
	// Compute diagonal.
	diagonal diag( compute_diagonal(p), closed_interval(p.y,p.y) );
	store_t::const_iterator i = store.lower_bound(diag);
	return i != store.end() && i->d == diag.d && i->r.contains(p.y);
}

int diagonal_store::compute_diagonal( point p ) const {
	// Slope is 1 and b = y - x. Thus When y == 0, x = -b.
	return -p.y + p.x;
}

void diagonal_store::clear() {
	store.clear();
}
