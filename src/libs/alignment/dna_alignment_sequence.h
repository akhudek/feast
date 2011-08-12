#ifndef DNA_ALIGNMENT_SEQUENCE_H_
#define DNA_ALIGNMENT_SEQUENCE_H_

#include <boost/shared_ptr.hpp>
#include <btl/sequence.h>
#include <iostream>
#include <btl/dna_gap.h>
#include <deque>

// standardize the name
typedef btl::dna_gap dna_alignment_alpha;

// dna sequence type
typedef std::deque<dna_alignment_alpha::symbol> dna_alignment_sequence_data;
typedef btl::sequence<dna_alignment_alpha,dna_alignment_sequence_data> dna_alignment_sequence;
typedef dna_alignment_sequence_data::iterator dna_alignment_sequence_data_iter;
typedef dna_alignment_sequence_data::const_iterator dna_alignment_sequence_data_citer;

// Smart Sequence Pointer
typedef boost::shared_ptr<dna_alignment_sequence> dna_alignment_sequence_ptr;

// allocator helper function
inline dna_alignment_sequence_ptr new_dna_alignment_sequence() {
   dna_alignment_sequence_ptr result = dna_alignment_sequence_ptr( new dna_alignment_sequence() );
   return result;
}  

inline std::ostream &write_dna_alignment_sequence_data( std::ostream &out, dna_alignment_sequence_data const &s ) {
	using namespace std;
	for( dna_alignment_sequence_data_citer i = s.begin(); i != s.end(); i++ ) {
		out << *i;
	}
	return out;
}

// vector of pointers
typedef std::vector<dna_alignment_sequence_ptr> dna_alignment_sequence_ptr_vector;
typedef dna_alignment_sequence_ptr_vector::iterator dna_alignment_sequence_ptr_vector_iter;

#endif /*DNA_SEQUENCE_H_*/
