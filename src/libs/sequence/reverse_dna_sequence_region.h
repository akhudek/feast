#ifndef REVERSE_DNA_SEQUENCE_REGION_H__
#define REVERSE_DNA_SEQUENCE_REGION_H__

#include "dna_sequence.h"
#include "reverse_sequence_region.h"

// dna sequence type
typedef reverse_sequence_region<dna_sequence_ptr> reverse_dna_sequence_region;
typedef reverse_dna_sequence_region::data_type reverse_dna_sequence_region_data;
typedef reverse_dna_sequence_region_data::iterator reverse_dna_sequence_region_data_iter;
typedef reverse_dna_sequence_region_data::const_iterator reverse_dna_sequence_region_data_citer;

// Smart Sequence Pointer
typedef boost::shared_ptr<reverse_dna_sequence_region> reverse_dna_sequence_region_ptr;

/*inline std::ostream &write_dna_sequence_region_data( std::ostream &out, dna_sequence_region_data const &s ) {
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
*/

#endif /*DNA_SEQUENCE_REGION_H_*/
