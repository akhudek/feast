#ifndef POINT_HASH_H__
#define POINT_HASH_H__
#include "point.h"

namespace __gnu_cxx {

	template<>
	struct hash<point> {
		size_t operator()( point const p ) const {
			int lowX = p.x & 0xFFFF;
			int lowY = p.y & 0xFFFF;
			return (size_t) (lowX << 16)| lowY;
		}
	};
}

#endif
