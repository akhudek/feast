#include "annotation.h"

std::ostream &operator<<( std::ostream &out, annotation const &a ) {
	for( annotation_citer i = a.begin(); i != a.end(); ++i ) {
		out << *i << std::endl;
	}
	return out;
}
