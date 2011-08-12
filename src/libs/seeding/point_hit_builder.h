/** Builds points from hits.
 *
 *  This assumes points are added to each sequence in the same order each time.
 *  Also, it assumes only two points are added to each hit.
 */

#ifndef POINT_HIT_BUILDER_H_
#define POINT_HIT_BUILDER_H_

#include "point.h"
#include "seed.h"

template<typename OutputIterator>
class point_hit_builder {
protected:
	enum builder_state { EMPTY= 0, FIRST_SET, SECOND_SET };
	// Internal state of builder.
	builder_state state;

	// Output Iterator
	OutputIterator oit;

	// Current point.
	point pnt;

public:
	point_hit_builder( OutputIterator out ) : state(EMPTY), oit(out), pnt(-1,-1) {}

	/** This actually does nothing in this case. */
	void start_hit();

	/** Add position in next sequence. */
	void add_position( int i, seed const &my_seed );

	/** Hit is finished, add it to oit. */
	void end_hit();
};

template<typename OutputIterator>
inline void point_hit_builder<OutputIterator>::start_hit() {}

template<typename OutputIterator>
inline void point_hit_builder<OutputIterator>::add_position( int i, seed const &my_seed ) {
	switch(state) {
	case EMPTY:
		pnt.x = i;
		state = FIRST_SET;
		break;
	case FIRST_SET:
		pnt.y = i;
		state = SECOND_SET;
		break;
	default:
		assert(0 && "point_hit_builder: adding to full hit");
	}
}

template<typename OutputIterator>
inline void point_hit_builder<OutputIterator>::end_hit() {
	assert( state == SECOND_SET );
	*oit = pnt;
	++oit;
	state = EMPTY;
}



#endif /* POINT_HIT_BUILDER_H_ */
