/**
 *  Structure for storing a 2d point.
 */
#ifndef POINT_H__
#define POINT_H__

#include <vector>
#include <iostream>
#include <cmath>

class point {
public:
	/** Location in first dimension. */
	int x;

	/** Location in second dimension. */
	int y;

	point() {}

	point(int new_a, int new_b ) : x(new_a), y(new_b) {}

	/** Test equality of two points. */
	bool operator==( point &o ) const;

	/** Rotates this point theta degrees around the origin. */
	void rotate( double theta );
};

inline void point::rotate( double theta ) {
	using namespace std;

	double oldx = (double)x, oldy = (double)y;
	x = (int)(cos(theta)*oldx - sin(theta)*oldy);
	y = (int)(sin(theta)*oldx + cos(theta)*oldy);
}

inline bool point::operator==( point &o ) const {
	return o.x == x && o.y == y;
}

// pretty output for closed interval
inline std::ostream &operator<<( std::ostream &out, point const &i ) {
	out << '(' << i.x << ',' << i.y << ')';
	return out;
}


#endif
