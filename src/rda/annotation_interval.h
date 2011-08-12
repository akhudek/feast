#ifndef ANNOTATION_INTERVAL_H__
#include "closed_interval.h"
#include <iostream>

class annotation_interval : public closed_interval {
	friend std::ostream &operator<<( std::ostream &out, annotation_interval const &ai );
public:
	int t;

	annotation_interval( int i, int j, int type );
};

#endif

