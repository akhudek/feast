#ifndef ALGORITHM_EXTRA_H__
#define ALGORITHM_EXTRA_H__
#include <algorithm>
#include <boost/random/uniform_int.hpp>
#include <boost/shared_ptr.hpp>

template<typename SequenceRegion>
void truncate_at_repeat( SequenceRegion &r ) {
	for( unsigned int i = 1; i < r.data.size(); ++i ) {
		if( r.data[i].masked() ) {
			r.truncate_at(i);
			break;
		}
	}
}

// filter the contents of the container by the given function
template< typename Container, typename FilterFunction >
void filter( Container &c, FilterFunction discard ) {
   typename Container::iterator new_end = std::remove_if( c.begin(), c.end(), discard );
   c.erase(new_end,c.end());
}

// apply function to each item
template< typename ForwardIterator, typename Function >
void apply( ForwardIterator begin, ForwardIterator end, Function f ) {
   for( ForwardIterator i = begin; i != end; i++ ) f(*i);
}

// like find, but returns the index of the position instead of an iterator
template< typename ForwardIterator, typename Value >
unsigned int index_of( ForwardIterator start, ForwardIterator end, Value const &v ) {
   unsigned int index = 0;
   for( ; start != end; start++, index++ ) if( *start == v ) return index;
   return index;
}

// split a container into two new containers
template< class Container, class FilterFunction, class OutputIteratorA, class OutputIteratorB >
void split_copy( Container &c, FilterFunction toA, OutputIteratorA outA, OutputIteratorB outB ) {
	for( typename Container::iterator i = c.begin(); i != c.end(); ++i ) {
		if( toA(*i) ) { *outA = *i; ++outA; }
		else { *outB = *i; ++outB; }
	}
}

// Uniformly sample entries.
template< typename Generator, class Container >
boost::shared_ptr<Container> uniform_sample( Generator &g, Container &c, unsigned int n ) {
	boost::shared_ptr<Container> newc( new Container() );
	for( unsigned int i = 0; i < n; ++i ) {
		boost::uniform_int<unsigned int> urng(0,c.size()-1);
		unsigned int j = urng( g );
		newc->push_back(c[j]);

		// swap out element
		c[j] = c.back();
		c.pop_back();
	}
	return newc;
}

/** Greater than or equal to according to second member of a pair.
 */
template< typename PairType >
struct gte_by2nd : public std::binary_function< PairType, PairType, bool > {
	gte_by2nd() {}

	inline bool operator()( PairType const &I,  PairType const &J) {
		return I.second >= J.second;
	}
};


#endif
