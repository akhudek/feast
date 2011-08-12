#ifndef CLOSED_INTERVAL_H__
#define CLOSED_INTERVAL_H__

#include <vector>
#include <iostream>

// Since we work with [start,end] ranges a lot, we define a type for them
// here.
class closed_interval {
public:
	int a;
	int b;

	closed_interval();

	closed_interval(int new_a, int new_b );

	/** Return the size of the interval. */
	int size() const;

	/** Returns true if point p is in this closed interval. */
	bool contains( int p ) const;

	bool operator==( closed_interval const &o ) const;
};

inline bool closed_interval::operator==( closed_interval const &o ) const {
	return a == o.a && b == o.b;
}


inline bool closed_interval::contains( int p ) const {
	return p >= a && p <= b;
}

inline int closed_interval::size() const {
	return b - a + 1;
}

typedef std::vector<closed_interval> closed_interval_vector;
typedef closed_interval_vector::iterator closed_interval_vector_iter;
typedef closed_interval_vector::const_iterator closed_interval_vector_citer;

// pretty output for closed interval
std::ostream &operator<<( std::ostream &out, closed_interval const &i );

#endif /*CLOSED_INTERVAL_H_*/
