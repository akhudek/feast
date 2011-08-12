#ifndef REVERSE_SEQUENCE_REGION_H_
#define REVERSE_SEQUENCE_REGION_H_
// Assumes Sequence uses a Random access container for storage.
#include <cassert>
#include "closed_interval.h"

template < typename Sequence_ptr >
class reverse_sequence_region {
private:
	// type of sequence
	typedef typename Sequence_ptr::element_type sequence;
	typedef typename sequence::data_type target_data_type;

public:
	typedef typename sequence::alphabet alphabet;

private:
	// This provides some of the vector interface, but really is just a pointer to
	// another vector.

	class vector_view {
	private:
		// the source sequence
		target_data_type *tdata;

	public:
		// the range within the source sequence
		closed_interval range;

		// iterators
		typedef std::reverse_iterator<typename target_data_type::const_iterator>	const_iterator;
		typedef std::reverse_iterator<typename target_data_type::iterator>	iterator;

		// empty data
		vector_view() {}

		// constructor
		vector_view( target_data_type &s, closed_interval r ) : tdata(&s), range(r) {}

		// copy constructor
		vector_view( vector_view const &o ) : tdata(o.tdata), range(o.range) {}

		// begining of the region
		inline const_iterator end() const {
			return const_iterator( tdata->begin() + range.a );
		}
		inline iterator end() {
			return iterator( tdata->begin() + range.a );
		}

		// end of the region
		inline const_iterator begin() const {
			return const_iterator( tdata->begin() + range.b + 1);
		}
		inline iterator begin() {
			return iterator( tdata->begin() + range.b + 1 );
		}

		inline typename alphabet::symbol &operator[]( unsigned int i ) {
			return (*tdata)[range.b-i];
		}

		// size of the region
		inline unsigned int size() {
			return range.size();
		}

		inline unsigned int size() const {
			return range.size();
		}

		vector_view &operator=( vector_view const &o ) {
			tdata = o.tdata;
			range = o.range;
			return *this;
		}

		/** Convert i to local coordinate for this region. */
		inline int local_coord( int i ) {
			return range.b - i;
		}

		/** Convert i to global coordinate for this region. */
		inline int global_coord( int i ) {
			return range.b - i;
		}

		inline closed_interval global_coord( closed_interval r ) {
			return closed_interval( range.b - r.b,  range.b - r.a );
		}

		inline closed_interval local_coord( closed_interval r ) {
			return closed_interval( range.b - r.b, range.b - r.a );
		}

	};

	inline closed_interval remap_range( closed_interval &a, closed_interval &b ) {
		closed_interval c;
		c.a = a.a + b.a;
		c.b = a.a + b.b;
		return c;
	}


public:
	// store a pointer to the original sequence to ensure it does not disappear.
	Sequence_ptr original_sequence;

	// supply an interface that is compatible with a sequence object
	typedef vector_view data_type;

	vector_view data;

	// Property map to hold sequence information.
	typedef typename sequence::tag_map tag_map;
	typedef typename tag_map::iterator tag_map_iter;
	typedef typename tag_map::const_iterator tag_map_const_iter;

	//tag_map tags;

	reverse_sequence_region() {}

	// constructor for the whole region
	reverse_sequence_region( Sequence_ptr s ) :
	original_sequence(s),
	data(s->data, closed_interval( 0, (int) s->data.size() - 1 ) ) {}
	//tags( s->tags ) {}

	// constructor
	reverse_sequence_region( Sequence_ptr s, closed_interval r ) : original_sequence(s), data( s->data, r ) { // , tags( s->tags ) {
		// check parameters
		assert( r.a >= 0 );
		assert( (unsigned int) r.b < s->data.size() );
	}

	// Constructor given an existing sequence region.
	reverse_sequence_region( reverse_sequence_region<Sequence_ptr> &s, closed_interval r ) : original_sequence(s.original_sequence),
	data( s.original_sequence->data, remap_range(s.data.range,r) ) {
	//tags( s.original_sequence->tags ) {
		// check parameters
		assert( r.a >= 0 );
		assert( (unsigned int) r.b < s.original_sequence->data.size() );
	}

	reverse_sequence_region &operator=( reverse_sequence_region const &o ) {
		data = o.data;
		//tags = o.tags;
		original_sequence = o.original_sequence;
		return *this;
	}

	reverse_sequence_region sub_region( closed_interval r ) {
		reverse_sequence_region subr( original_sequence, closed_interval( data.range.b - r.b, data.range.b - r.a ) );
		return subr;
	}

	void truncate_at( int i ) {
		assert( i > 0 );
		assert( i < data.range.size() );
		data.range.a = data.range.b - i + 1;
	}

};

#endif
