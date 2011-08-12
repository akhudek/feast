#include "extension_store.h"
#include <boost/tuple/tuple.hpp>

extension_store::extension_store() {}

extension_store::extension::row::row() : range(-1,-1), best_col(-1) {}
extension_store::extension::row::row( closed_interval r, int c ) : range(r), best_col(c) {}

extension_store::extension::extension() : range(-1,-1) {}

void extension_store::clear() {
	store.clear();
}

extension_store::extension_ptr extension_store::split( extension_store::extension_ptr e, int j ) {
	using namespace std;
	assert( j >= e->range.a );

	// Create tail and link to this extension.
	int offset = j - e->range.a;

	extension_ptr T( new extension() );
	copy( e->fill.begin() + offset, e->fill.end(), back_inserter( T->fill ) );
	T->range.a = j;
	T->range.b = e->range.b;

	// Truncate this extension.
	e->fill.resize( offset );
	e->range.b = j - 1;

	return T;
}

extension_store::store_map::iterator extension_store::truncate_split( extension_store::store_map::iterator t, int j ) {
	using namespace std;
	// Since segments are stored by the lower range, we need to leave the lower part of the segment.

	// We insert a new segment for the truncated portion.
	store_map::iterator new_t = store.insert( t, make_pair( j-1, segment_map() ) );
	segment_map::iterator new_s = new_t->second.end();
	int offset = j - t->second.begin()->second->range.a;

	for( segment_map::iterator s = t->second.begin(); s != t->second.end(); ++s ) {
		// Create new extension s and copy the first values into it.
		extension_ptr new_e( new extension() );
		new_e->range.a = s->second->range.a;
		new_e->range.b = j-1;
		copy( s->second->fill.begin(), s->second->fill.begin() + offset, back_inserter( new_e->fill ) );

		// Adjust values in old segment.
		s->second->range.a = j;
		copy( s->second->fill.begin()+offset, s->second->fill.end(), s->second->fill.begin() );
		s->second->fill.resize( s->second->range.size() );

		// Add new segment.
		new_s = new_t->second.insert( new_s, make_pair( new_e->fill.back().range.b , new_e) );

		// Adjust next and prev.
		new_e->prev = s->second->prev;
		// Update next of prev if it exists.
		if( new_e->prev ) new_e->prev.get()->second->next = new_s;
		new_e->next = s;
		s->second->prev = new_s;
	}
	return new_t;
}

boost::optional<extension_store::segment_map::iterator> extension_store::add( extension_store::extension_ptr e ) {
	using namespace std;
	using namespace boost;

	// Find first segment such that range.b is a lower_bound on range.a.
	store_map::iterator t = store.lower_bound( e->range.a );
	if( t == store.end() ) return insert( t, e );

	segment_map::iterator s = t->second.begin();
	if( s->second->range.a > e->range.b ) return insert( t, e );

	// Otherwise we conflict with with one or more segments. We have two
	// possible cases.

	// Here we need to split e at s.range.a.
	//  [     e    ]        ->   [ e ][  new_e  ]
	//     [      s     ]	          [     s          ]
	if( s->second->range.a > e->range.a ) {
		extension_ptr new_e = split(e, s->second->range.a);
		new_e->prev = insert( t, e );	// e does not conflict with anything, so add it.
		e->next = add( new_e );			// Deal with the remaining segment and update next.
		return new_e->prev;
	}

	// Here we need to split s and all it's friends at e.range.a.
	//       [     e    ]        ->        [    e    ]
	//  [      s              ]	     [    ][    s         ]
	//							or
	//       [    e     ]		 ->        [    e    ]
	//  [     s    ]		         [    ][  s  ]
	//							or
	//       [    e     ]		 ->        [    e    ]
	//  [     s         ]            [    ][  s      ]

	if( s->second->range.a < e->range.a ) {
		truncate_split( t, e->range.a );
	}

	// Now we have s.range.a == e.range.a with s pointing to the leftmost segment
	// in a set. We have three possible cases.

	// We have s shorter than e, so we split e at (**s).range.b+1
	//  [     e    ]   ->   [  e  ][  new_e  ]
	//  [  s  ]	            [  s  ]
	if( s->second->range.b < e->range.b ) {
		extension_ptr new_e = split(e, s->second->range.b+1);
		e->next = add( new_e );
	}

	// We have s longer than e, so we need to split s again.
	//  [     e     ]          ->    [     e     ]
	//  [     s            ]         [     s     ][    ]
	if( s->second->range.b > e->range.b ) {
		t = truncate_split( t, e->range.b + 1 );
		s = t->second.begin();
	}

	//  [    e     ]
	//  [    s     ]
	assert( s->second->range == e->range  );			// Now we must have s.range == e.range.
	segment_map::iterator new_s = insert( t, e );		// We insert e into t.
	if( e->next ) e->next.get()->second->prev = new_s;  // Update prev in next if it exists.
	return new_s;
}

extension_store::bounds extension_store::get_bounds( point p, bool d ) {
	using namespace std;
	using namespace boost;

	bounds b( store, d );

	b.j = p.y;
	b.t = store.lower_bound(p.y);
	b.rt = store_map::const_reverse_iterator(b.t);

	// No possible segment with match.
	if( b.t == store.end() ) {
		// I have no idea why this only works with b.store. It must be because of const.
		if( !b.forward && b.rt != b.store.rend() ) {
			segment_map::const_iterator s = b.rt->second.begin();
			b.j_begin = s->second->range.b;
		}
		return b;
	// Reverse iterators are always one previous of the iterator they are constructed on.
	// So we move it back to match t.
	} else --b.rt;

	segment_map::const_iterator s = b.t->second.begin();

	// Segment s overlaps j.
	if( s->second->range.a <= p.y ) {
		b.s_j = p.y - s->second->range.a;
		b.s_end = s->second->fill.size();

		segment_map::const_iterator r = b.t->second.upper_bound( p.x );
		tie(b.L,b.R) = search_around(r,b.t->second,b.s_j,closed_interval(p.x,p.x) );

		// Check for point covered.
		if( (b.L && b.L->fill[b.s_j].range.b >= p.x ) || (b.R && b.R->fill[b.s_j].range.a <= p.x) )
			throw point_covered_exception();

		return b;
	}

	// Otherwise segment does not overlap j, so we continue.
	if( b.forward ) b.j_begin = s->second->range.a;
	else if( b.rt != b.store.rend() ) {
		++b.rt;
		if( b.rt != b.store.rend() ) b.j_begin = b.rt->second.begin()->second->range.b;
		}

	//cerr << b.t->second.begin()->second->range << "\t" << b.j << "\t" << b.j_begin << "\t" << b.s_j << "\t" << b.s_end << endl;

	return b;
}


extension_store::bounds::bounds( store_map const &s, bool d )
	: store(s), forward(d),  s_j(-1), j_begin(-1), s_end(-1), j(-1) {}


void extension_store::check() {
	using namespace std;
	// First ensure that all segments in each y-range have the correct y-range.
	int last = -1;
	for( extension_store::store_map::iterator t = store.begin(); t != store.end(); ++t ) {
		assert( t->second.size() > 0 );
		cerr << "T " << t->second.begin()->second->range << "\t" << t->first << endl;

		for( segment_map::iterator s = t->second.begin(); s != t->second.end(); ++s ) {
			cerr << "S " << s->second->range << endl;

			assert( s->second->fill.size() == (unsigned int)s->second->range.size() );
			assert( s->second->range.b == t->first );
			assert( s->second->range.b == t->second.begin()->second->range.b );
			assert( s->second->range.a == t->second.begin()->second->range.a );
			assert( s->second->range.a > last );

			segment_map::iterator s2 = s;
			++s2;
			if( s2 != t->second.end() ) {
				for( unsigned int i = 0; i < s->second->fill.size(); ++i ) {
					assert( s->second->fill[i].range.a <= s->second->fill[i].range.b );
					assert( s->second->fill[i].range.b < s2->second->fill[i].range.a );
				}
			}


			if( s->second->next ) assert( s->second->next.get()->second->prev.get() == s );
			if( s->second->prev ) assert( s->second->prev.get()->second->next.get() == s );
		}

		last = t->first;
	}

}





