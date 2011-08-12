#include "annotation_interval.h"

annotation_interval::annotation_interval( int i, int j, int k ) : t(k) {
	a = i;
	b = j;
}

std::ostream &operator<<( std::ostream &out, annotation_interval const &i ) {
	out << "[ " << i.a << " " << i.b << " ] " << i.t;
	return out;
}
