#include "chaining_algorithm.h"

std::vector<local_chain>::iterator get_chain_middle_it( std::vector<local_chain>::iterator i, std::vector<local_chain>::iterator end ) {
	using namespace std;
	int total_pnt = 0, mcnt = 1;
	vector<local_chain>::iterator j = i;
	while( i != end ) {
		// Move forward in i.
		++total_pnt;
		i = i-> previous;

		// If we need to advance middle, do so.
		while( mcnt < total_pnt/2 ) {
			++mcnt;
			j = j->previous;
		}
	}
	return j;
}

point get_chain_middle( std::vector<local_chain>::iterator i, std::vector<local_chain>::iterator end ) {
	return get_chain_middle_it(i,end)->p;
}
