#ifndef FRAGMENT_H__
#define FRAGMENT_H__

#include <iostream>
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>
#include "closed_interval.h"

// A fragment is a set of ranges related by some criteria.
class fragment {
   // allow functions to access private data
   friend bool operator<<( fragment const &a, fragment const &b );
   friend bool operator>>( fragment const &a, fragment const &b );
   friend std::ostream &operator<<( std::ostream &out, fragment const &f );

protected:
	std::vector<closed_interval> regions;

public:
   fragment() {}

	fragment( fragment const &o ) {
		regions = o.regions;
	}


   // add an interval
   void add( closed_interval c );

   // number of ranges
   unsigned int size() const;

   // return the region at position i
   closed_interval const &operator[]( unsigned int i ) const;
   closed_interval &operator[]( unsigned int i );

};

inline closed_interval &fragment::operator[]( unsigned int i )  {
   return regions[i];
}


inline closed_interval const &fragment::operator[]( unsigned int i ) const {
   return regions[i];
}

inline void fragment::add( closed_interval c ) {
   regions.push_back( c );
}

inline unsigned int fragment::size() const {
   return regions.size();
}

// define some shared pointers
typedef boost::shared_ptr<fragment> fragment_ptr;

// define fragment vectors
typedef std::vector<fragment_ptr> fragment_ptr_vector;
typedef fragment_ptr_vector::iterator fragment_ptr_vector_iter;
typedef fragment_ptr_vector::reverse_iterator fragment_ptr_vector_riter;
typedef fragment_ptr_vector::const_iterator fragment_ptr_vector_citer;
typedef boost::shared_ptr<fragment_ptr_vector> fragment_ptr_vector_ptr;

// define fragment lists
typedef std::list<fragment_ptr> fragment_ptr_list;
typedef fragment_ptr_list::iterator fragment_ptr_list_iter;
typedef fragment_ptr_list::const_iterator fragment_ptr_list_citer;
typedef boost::shared_ptr<fragment_ptr_list> fragment_ptr_list_ptr;

#endif
