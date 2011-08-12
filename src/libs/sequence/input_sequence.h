#ifndef INPUT_SEQUENCE_H_
#define INPUT_SEQUENCE_H_

#include <boost/shared_ptr.hpp>
#include <btl/sequence.h>
#include <btl/mdna_iupac.h>
#include <iostream>

// we allow iupac input
typedef btl::mdna_iupac input_alpha;

// input sequence type
typedef btl::sequence<input_alpha> input_sequence;
typedef input_sequence::data_type input_sequence_data;
typedef input_sequence_data::iterator input_sequence_data_iter;
typedef input_sequence_data::const_iterator input_sequence_data_citer;

// Smart Sequence Pointer
typedef boost::shared_ptr<input_sequence> input_sequence_ptr;

inline std::ostream &write_input_sequence_data( std::ostream &out, input_sequence_data const &s ) {
	using namespace std;
	copy(s.begin(),s.end(),ostream_iterator<input_alpha::symbol>(out));
	return out;
}
  
// vector of pointers
typedef std::vector<input_sequence_ptr> input_sequence_ptr_vector;
typedef input_sequence_ptr_vector::iterator input_sequence_ptr_vector_iter;

#endif /*INPUT_SEQUENCE_TYPES_H_*/
