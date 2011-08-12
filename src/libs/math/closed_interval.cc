#include "closed_interval.h"

closed_interval::closed_interval() : a(0), b(0) {}
		
closed_interval::closed_interval(int new_a, int new_b ) : a(new_a), b(new_b) {}

// pretty output for closed interval
std::ostream &operator<<( std::ostream &out, closed_interval const &i ) {
	out << "[ " << i.a << ' ' << i.b << " ]";
	return out;
}
