/**
 * Hashes a sequence into a hash data class according to a seed. This uses the new
 * fast_hash_hit_data interface.
 */

#ifndef FAST_HASH_FACTORY_H_
#define FAST_HASH_FACTORY_H_

#include <iostream>
#include "seed.h"
#include "sequence_region.h"

template<class HashResult>
class fast_hash_factory {
protected:
	/** Seed we hash with. */
	 seed my_seed;

	/** Use masking information. */
	bool use_mask;

	/** Hash a sequence into the hash data using seed s. */
	template< typename SequenceRegion >
	void hash_sequence( HashResult &data, SequenceRegion &seq );

	/** Hash data starting at i to data according to my_seed.
	 * Here i is assumed to be at position pos.
	 **/
	template< typename ConstIterator >
	ConstIterator hash_position( HashResult &data, ConstIterator i, ConstIterator end, unsigned int pos );

public:

	//--------------------------------------------------------------------------
	// Constructors.
	//--------------------------------------------------------------------------
	fast_hash_factory( seed &new_seed, bool mask = false );

	/** Hash each sequence region from begin to end.
	 *  We assume each iterator dereferences to a pointer to a sequence region of
	 *  the correct type.
	 */
	template< typename ConstIterator >
	boost::shared_ptr<HashResult> create_from_ptr( ConstIterator begin, ConstIterator end );

};


template<class HashResult>
fast_hash_factory<HashResult>::fast_hash_factory(seed &new_seed, bool mask) :
	my_seed(new_seed), use_mask(mask) {}

template<class HashResult>
 template<typename ConstIterator>
boost::shared_ptr<HashResult> fast_hash_factory<HashResult>::create_from_ptr( ConstIterator it_start, ConstIterator it_end ) {
	// create new hash and set match length
	boost::shared_ptr<HashResult> data( new HashResult(my_seed) );

	// hash each sequence
	for( ConstIterator i = it_start; i != it_end; ++i) hash_sequence( *data, **i );
	return data;
}


template<class HashResult>
 template<typename SequenceRegion>
void fast_hash_factory<HashResult>::hash_sequence( HashResult &data, SequenceRegion &seq ) {
	using namespace std;

	data.start_new_sequence();

	// break if region is too small to use seed
	if( seq.data.size() < my_seed.length() ) return;

	// find last possible match position + 1 (prevent seed from going past end)
	typename SequenceRegion::data_type::const_iterator end, real_end = seq.data.end();
	end = seq.data.begin() + seq.data.size() - my_seed.length() + 1;

	// Use position information from region.
	unsigned int pos = seq.data.range.a;

	for( typename SequenceRegion::data_type::const_iterator i = seq.data.begin(), j; i != end; ) {
		// hash position i
		j = hash_position( data, i, real_end, pos );
		if( j == i ) {
			++pos;
			++i;
		// skip masked area
		} else {
			pos += j - i;
			i = j;
			// Check for the case where by skipping masked sequence we end up in  [end,real_end).
			if( real_end - i < (int) my_seed.length() ) break;
		}
	}
}

template<class HashResult>
 template<typename ConstIterator>
ConstIterator fast_hash_factory<HashResult>::hash_position( HashResult &data, ConstIterator it, ConstIterator end, unsigned int p ) {
	ConstIterator i = it;
	typename HashResult::hash_key k = 0;
	for( seed::seed_mask::const_iterator j = my_seed.pattern_mask().begin(); j != my_seed.pattern_mask().end(); ++j, ++i ) {
		// Encountered masked base. Skip to next good base and do not hash.
		if( i->masked() ) {
			while( i != end && i->masked() ) ++i;
			return i;
		// Otherwise, if we are in a 1 position add symbol to key.
		} else if( *j ) {
			// We assume that a DNA base is encoded in two bits.
			k <<= 2;
			k |= i->index();
		}
	}
	data.add( k, p );
	return it;
}

#endif /* FAST_HASH_FACTORY_H_ */
