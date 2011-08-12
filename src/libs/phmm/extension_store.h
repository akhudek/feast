/** Store alignment extensions for two sequences.
 *
 *  This data structure is organized as follows. The class extension stores
 *  a part of an alignment extension covering a range of rows. It stores a
 *  pointer called next to the next part of the extension, if needed.
 *
 *  The binary tree store, stores extension segments ordered by range.b. All
 *  extensions in store have the following properties.
 *
 *  1) If extension segments A and B overlap, they overlap completely. That
 *     is A.range == B.range.
 *
 *  2) Overlapping extension segments never cross. For extension segments
 *     A and B where A.range == B.range and A.col < B.col, we always have
 *     A.fill[i].region.b < B.fill[i].region.a for all i.
 *
 *  3) Every set of overlapping extension segments A1, ..., An, is linked
 *     by left and right pointers and ordered by Ai.col. The segment
 *     pointed to by store is always the left most segment.
 */

#ifndef EXTENSION_STORE_H_
#define EXTENSION_STORE_H_

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <map>
#include <btl/logspace.h>
#include "closed_interval.h"
#include "point.h"

#include "ereal.h"

class extension_store {
public:
	/** Data structure for storing extension segments. */
	struct extension {
		struct row {
			closed_interval range;
			int best_col;
			row();
			row( closed_interval r, int c );
		};
		std::vector<row> fill;
		closed_interval range;

		boost::optional< std::map< int, boost::shared_ptr<extension> >::iterator > next;
		boost::optional< std::map< int, boost::shared_ptr<extension> >::iterator > prev;

		void add_row( closed_interval r, int c );

		extension();
	};
	typedef boost::shared_ptr<extension> extension_ptr;

	typedef std::map<int,extension_ptr> segment_map;
	typedef std::map<int,segment_map>	store_map;

protected:
	store_map store;

	segment_map::iterator insert( store_map::iterator t, extension_ptr e );

	/** Split extension into two new extensions.
	 *  Given extension over 'range' and split point j, we return two new
	 *  extensions A and B such that A.range = (range.a,j-1) and
	 *  B.range = (j,range.b).
	 */
	extension_ptr split( extension_ptr e, int j );
	store_map::iterator truncate_split( store_map::iterator t, int j );


public:
	class bounds {
	protected:
		friend class extension_store;

		store_map const &store;

		bool forward;

		/** Current or next segment. */
		extension_store::store_map::const_iterator t;
		extension_store::store_map::const_reverse_iterator rt;

		/** Current position in segment. */
		int s_j;

		int j_begin;

		/** End of current segment. */
		int s_end;

		extension_ptr L;

		extension_ptr R;

		/** Current position overall. */
		int  j;


	public:
		bounds( extension_store::store_map const &s, bool d );
		closed_interval clip( closed_interval b );
		void step( closed_interval b );

		bounds &operator=( bounds const &o );
	};

	/** Add a new extension to the store, we assume it satisfies
	 *  property (2) after split into appropriate segments.
	 */
	boost::optional<segment_map::iterator> add( extension_ptr );

	bounds get_bounds( point x, bool d );

	/** Search around point. We have several cases.
	 *
	 *  A)   [     r     ]
	 *  					[     p     ]
	 *
	 *
	 *  B)   [     r     ]
	 *   			[     p     ]
	 *
	 *
	 *  C)   [     r     ]
	 *  	   [   p   ]
	 *
	 *
	 *  D)        [     r     ]
	 *       [         p         ]
	 *
	 *
	 *  E)        [     r     ]
	 *       [     p     ]
	 *
	 *
	 *  F)                    [     r     ]
	 *       [     p     ]
	 *
	 *  Both B and C mean we should probably just give up. There are unlikely to
	 *  be any extreme cases that are interesting. We arbitrarily assign these
	 *  as left segments and scan right.
	 */

	static std::pair<extension_ptr,extension_ptr> search_around( segment_map::const_iterator s, segment_map const &t, int s_j, closed_interval r );
	static std::pair<extension_ptr,extension_ptr> search_left( segment_map::const_iterator s, segment_map const &t, int s_j, closed_interval r );
	static std::pair<extension_ptr,extension_ptr> search_right( segment_map::const_iterator s, segment_map const &t, int s_j, closed_interval r );

	extension_store();

	void clear();
	void check();
};


inline void extension_store::extension::add_row( closed_interval r, int c ) {
	fill.push_back( row( r, c) );
}

inline extension_store::segment_map::iterator extension_store::insert( extension_store::store_map::iterator t, extension_ptr e ) {
	using namespace std;
	store_map::iterator new_t = store.insert( t, make_pair( e->range.b, segment_map() ) );
	return (new_t->second.insert( make_pair( e->fill.back().range.b, e ) )).first;
}


inline std::pair<extension_store::extension_ptr,extension_store::extension_ptr>
	extension_store::search_around( extension_store::segment_map::const_iterator s, extension_store::segment_map const &t, int s_j, closed_interval r ) {
	assert( t.size() > 0 );

	if( s == t.end() ) --s;

	closed_interval &p = s->second->fill[s_j].range;

	if( p.b > r.b && p.a >= r.a ) return search_left(s,t,s_j,r);

	return search_right(s,t,s_j,r);
}

inline std::pair<extension_store::extension_ptr,extension_store::extension_ptr>
	extension_store::search_left( extension_store::segment_map::const_iterator s, extension_store::segment_map const &t, int s_j, closed_interval r ) {
	using namespace std;
	extension_ptr L, R;
	while(1) {

		closed_interval &p = s->second->fill[s_j].range;

		// Case A or B.
		if( p.b > r.b && p.a >= r.a ) R = s->second;
		else L = s->second;

		if( p.b < r.a || s == t.begin() ) break;

		--s;
	}
	return make_pair(L,R);
}

inline std::pair<extension_store::extension_ptr,extension_store::extension_ptr>
	extension_store::search_right( extension_store::segment_map::const_iterator s, extension_store::segment_map const &t, int s_j, closed_interval r ) {
	using namespace std;
	extension_ptr L, R;

	while(1) {

		closed_interval &p = s->second->fill[s_j].range;

		// Case A or B.
		if( p.b > r.b && p.a >= r.a ) R = s->second;
		else L = s->second;

		++s;

		if( p.a > r.b || s == t.end() ) break;

	}
	return make_pair(L,R);
}




class point_covered_exception : public std::exception {
	virtual const char* what() const throw() {
		return "Point covered by existing extension.";
	}
};

inline extension_store::bounds &extension_store::bounds::operator=( extension_store::bounds const &o ) {
	// We assume store is the same.
	assert( &store == &o.store );
	forward = o.forward;
	t = o.t;
	rt = o.rt;
	s_j = o.s_j;
	j_begin = o.j_begin;
	s_end = o.s_end;
	L = o.L;
	R = o.R;
	j = o.j;
	return *this;
}

inline void extension_store::bounds::step( closed_interval next_b ) {
	using namespace boost;
	using namespace std;
	if( forward ){
		if( t == store.end() ) return;	// We'll never have anything to do.
		++j;
		if( s_j > -1 ) {				// We are still in the segment.
			++s_j;
			if( s_j == s_end ) {		// We reached the end.
				s_j = -1;

				++t;					// Move t up.

				if( t == store.end() ) return;

				// If we have a next L segment, use to set L and R.
				if( L && L->next ) {

					s_j = 0;
					s_end = t->second.begin()->second->fill.size();
					R.reset();
					segment_map::const_iterator r = *(L->next);
					tie(L,R) = extension_store::search_right(r,t->second,s_j,next_b);

					// If we have a next R segment, use to set L and R.
				} else if( R && R->next ) {

					s_j = 0;
					s_end = t->second.begin()->second->fill.size();
					L.reset();
					segment_map::iterator r = *(R->next);
					tie(L,R) = extension_store::search_left(r,t->second,s_j,next_b);

				// No L or R, but we may still have a new segment.
				// We leave this case for the code below.
				} else {
					L.reset();
					R.reset();
					j_begin = t->second.begin()->second->range.a;
				}
			}
		}

		// We need to start an alignment.
		if( s_j == -1 && j == j_begin ) {

			s_j = 0;
			s_end = t->second.begin()->second->fill.size();
			segment_map::const_iterator r = t->second.upper_bound( next_b.a );
			tie(L,R) = extension_store::search_around(r,t->second,s_j,next_b );
		}
		// Otherwise, we are waiting for the next segment.

	// Otherwise we are moving backwards.
	} else {
		if( rt == store.rend() ) return;	// We'll never have anything to do.
		--j;
		if( s_j > -1 ) {				// We are still in the segment.
			--s_j;
			if( s_j == -1 ) {			// We reached the end.

				++rt;					// Move t up.

				if( rt == store.rend()  ) return;

				// If we have a prev L segment, use to set L and R.
				if( L && L->prev ) {
					s_j = rt->second.begin()->second->fill.size()-1;
					R.reset();
					segment_map::iterator r = *(L->prev);

					tie(L,R) = extension_store::search_right(r,rt->second,s_j,next_b);

				// If we have a prev R segment, use to set L and R.
				} else if( R && R->prev ) {
					assert( rt-> second.size() > 0 );
					s_j = rt->second.begin()->second->fill.size() - 1;
					L.reset();
					segment_map::iterator r = *(R->prev);

					tie(L,R) = extension_store::search_left(r,rt->second,s_j,next_b);

				// No L or R, but we may still have a new segment.
				// We leave this case for the code below.
				} else {
					L.reset();
					R.reset();
					j_begin = rt->second.begin()->second->range.b;
				}
			}
		}

		// We need to start an alignment.
		if( s_j == -1 && j == j_begin ) {
			s_j = rt->second.begin()->second->fill.size()-1;
			segment_map::const_iterator r = rt->second.upper_bound( next_b.a );
			tie(L,R) = search_around(r,rt->second,s_j,next_b );
		}
		// Otherwise, we are waiting for the next segment.
	}
}

inline closed_interval extension_store::bounds::clip( closed_interval r ) {
	using namespace std;
	// Nothing to clip.
	if( s_j == -1 || (forward && t == store.end()) || (!forward && rt == store.rend()) ) return r;

	if( R ) r.b = min( r.b, R->fill[s_j].range.a - 1 );
	if( L ) r.a = max( r.a, L->fill[s_j].range.b + 1 );

	return r;
}


#endif /* EXTENSION_STORE_H_ */
