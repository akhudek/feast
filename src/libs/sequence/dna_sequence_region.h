#ifndef DNA_SEQUENCE_REGION_H__
#define DNA_SEQUENCE_REGION_H__

#include "dna_sequence.h"
#include "sequence_region.h"

// dna sequence type
typedef sequence_region<dna_sequence_ptr> dna_sequence_region;
typedef dna_sequence_region::data_type dna_sequence_region_data;
typedef dna_sequence_region_data::iterator dna_sequence_region_data_iter;
typedef dna_sequence_region_data::const_iterator dna_sequence_region_data_citer;

// Smart Sequence Pointer
typedef boost::shared_ptr<dna_sequence_region> dna_sequence_region_ptr;

inline std::ostream &write_dna_sequence_region_data( std::ostream &out, dna_sequence_region_data const &s ) {
	using namespace std;
	for( dna_sequence_region_data_citer i = s.begin(); i != s.end(); i++ ) {
		out << *i;
	}
	return out;
}

// allocator
inline dna_sequence_region_ptr new_dna_sequence_region( dna_sequence_ptr s ) {
   dna_sequence_region_ptr newptr( new dna_sequence_region( s ) );
   return newptr;
}

// vector of pointers
typedef std::vector<dna_sequence_region_ptr> dna_sequence_region_ptr_vector;
typedef dna_sequence_region_ptr_vector::iterator dna_sequence_region_ptr_vector_iter;

#endif /*DNA_SEQUENCE_REGION_H_*/
