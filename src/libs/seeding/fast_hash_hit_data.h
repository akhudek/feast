/** Faster but less memory efficient version of hash_hit_data.
 *
 *	In each bin we store a singly linked list of positions and a
 *	list of pointers into the list. The pointers indicate the
 *	end of positions from one sequence and start of another.
 */

#ifndef FAST_HASH_HIT_DATA_H_
#define FAST_HASH_HIT_DATA_H_

#include <vector>
#include <ext/slist>
#include <ext/hash_map>
#include "null_hash.h"
#include "seed.h"

template< typename ALPHABET >
class fast_hash_hit_data {
public:
	/** Data type for position list. */
	typedef __gnu_cxx::slist<unsigned int> position_list;

	/** Data type for sequence ends in position list. */
	typedef __gnu_cxx::slist<position_list::const_iterator> sequence_end_list;

	struct position_data {
		sequence_end_list sequence_ends;
		position_list positions;

	};
	typedef unsigned int hash_key;
	typedef __gnu_cxx ::hash_map< hash_key, position_data, null_hash<hash_key> > hash_data;

protected:
	hash_data hash;
	unsigned int num_sequences;
	seed 	  my_seed;

public:
	fast_hash_hit_data( seed &newseed );

	/** Return a reference to the hashed data. */
	hash_data const &get_hash() const;

	/** Add pos to hash with key. */
	void add( hash_key const &key, unsigned int pos );

	/** Add a new end of sequence marker to every bin. */
	void start_new_sequence();

	/** Return number of sequences. */
	unsigned int number_of_sequences() const;

	/** Return seed hash was created with. */
	seed const &get_seed() const;

	/** Enumerate hits into objects.
	 *  What objects are create depend on the HitBuilder class.
	 */
	template< class HitBuilder >
	void enumerate( HitBuilder hbuild );
};

// pointer definition

template< typename ALPHABET >
fast_hash_hit_data<ALPHABET>::fast_hash_hit_data( seed &newseed ) : num_sequences(0), my_seed(newseed) {}

template< typename ALPHABET >
inline typename fast_hash_hit_data<ALPHABET>::hash_data const &fast_hash_hit_data<ALPHABET>::get_hash() const {
	return hash;
}

template< typename ALPHABET >
inline seed const &fast_hash_hit_data<ALPHABET>::get_seed() const {
	return my_seed;
}

template< typename ALPHABET >
inline unsigned int fast_hash_hit_data<ALPHABET>::number_of_sequences() const {
	return num_sequences;
}

template< typename ALPHABET >
inline void fast_hash_hit_data<ALPHABET>::add( typename fast_hash_hit_data<ALPHABET>::hash_key const &key, unsigned int pos ) {
	position_data &hash_bin = hash[key];
	// If we haven't seen this bin before, we must add the correct number end pointers.
	if( hash_bin.sequence_ends.empty() )
		for( unsigned int i = 0; i < num_sequences; ++i )
			hash_bin.sequence_ends.push_front( hash_bin.positions.end() );
	hash_bin.positions.push_front(pos);
}

template< typename ALPHABET >
inline void fast_hash_hit_data<ALPHABET>::start_new_sequence() {
	// First element is now end of new sequence.
	for( typename hash_data::iterator i = hash.begin(); i != hash.end(); ++i )
		i->second.sequence_ends.push_front( i->second.positions.begin() );

	++num_sequences;
}

template< typename ALPHABET >
template< class HitBuilder >
inline void fast_hash_hit_data<ALPHABET>::enumerate( HitBuilder builder ) {
	using namespace std;

	vector<position_list::const_iterator> position_breaks( num_sequences + 1 ),	cur_position( num_sequences );

	for( typename hash_data::const_iterator i = hash.begin(); i != hash.end(); i++) {
		// ensure hits in all sequences and save start positions
		bool complete = true;

		// if last sequence is empty then end matches beginning
		if( *(i->second.sequence_ends.begin()) == i->second.positions.begin() ) continue;

		// setup up pointers to sequence break points
		vector<position_list::const_iterator>::iterator pb = position_breaks.begin(), cp;
		*pb = i->second.positions.begin();
		++pb;
		*pb = i->second.sequence_ends.front();
		++pb;

		// successive end points indicate empty position list
		for( sequence_end_list::const_iterator j = i->second.sequence_ends.begin(),
			 j_next = ++(i->second.sequence_ends.begin());
			 j_next != i->second.sequence_ends.end(); ++j, ++j_next, ++pb ) {

			if( *j == *j_next ) complete = false;
			*pb = *j_next;
		}
		if( !complete ) continue;

		// complete position, now we enumerate
		for(cp = cur_position.begin(), pb = position_breaks.begin(); cp != cur_position.end(); ++cp, ++pb ) *cp = *pb;
		vector<position_list::const_iterator>::reverse_iterator rpb, rcp;

		while( cur_position.front() != *(++position_breaks.begin()) ) {

			builder.start_hit();

			for(rcp = cur_position.rbegin(); rcp != cur_position.rend(); ++rcp ) {
				builder.add_position( **rcp, my_seed );
			}

			builder.end_hit();

			// move cp pointers
			rpb = position_breaks.rbegin(), rcp = cur_position.rbegin();
			++(*rcp);
			while( *rcp == *rpb && cur_position.front() != *(++position_breaks.begin()) ) {
				++rpb;
				*rcp = *rpb;
				++rcp;
				++(*rcp);
			}
		}
	}
}


#endif /* FAST_HASH_HIT_DATA_H_ */
