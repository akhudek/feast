#ifndef DNA_SEQUENCE_H_
#define DNA_SEQUENCE_H_

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <btl/sequence.h>
#include <btl/mdna.h>

// standardize the name
typedef btl::mdna dna_alpha;

// dna sequence type
typedef btl::sequence<dna_alpha> dna_sequence;
typedef dna_sequence::data_type dna_sequence_data;
typedef dna_sequence_data::iterator dna_sequence_data_iter;
typedef dna_sequence_data::const_iterator dna_sequence_data_citer;

// Smart Sequence Pointer
typedef boost::shared_ptr<dna_sequence> dna_sequence_ptr;

// allocator helper function
inline dna_sequence_ptr new_dna_sequence() {
   dna_sequence_ptr result( new dna_sequence() );
   return result;
}  

inline std::ostream &write_dna_sequence_data( std::ostream &out, dna_sequence_data const &s ) {
	using namespace std;
	for( dna_sequence_data_citer i = s.begin(); i != s.end(); i++ ) {
		out << *i;
	}
	return out;
}

// vector of pointers
typedef std::vector<dna_sequence_ptr> dna_sequence_ptr_vector;
typedef dna_sequence_ptr_vector::iterator dna_sequence_ptr_vector_iter;
typedef dna_sequence_ptr_vector::const_iterator dna_sequence_ptr_vector_citer;

// smart pointer for vector of pointers
typedef boost::shared_ptr<dna_sequence_ptr_vector> dna_sequence_ptr_vector_ptr;

// allocator helper function
inline dna_sequence_ptr_vector_ptr new_dna_sequence_ptr_vector() {
   dna_sequence_ptr_vector_ptr result( new dna_sequence_ptr_vector() );
   return result;
}  


#endif /*DNA_SEQUENCE_H_*/
