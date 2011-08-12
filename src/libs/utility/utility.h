/*
 * utility.h
 *
 *  Created on: 23-Jun-2009
 *      Author: akhudek
 */

#ifndef UTILITY_H_
#define UTILITY_H_
#include "fragment.h"
#include "dna_sequence_region.h"
#include "io.h"

template<class SEQ>
void unmask( SEQ &s ) {
   for( typename SEQ::data_type::iterator i = s.data.begin(); i != s.data.end(); ++i ) i->set_masked(false);
}

template<class SEQ>
bool is_all_masked( SEQ &s ) {
   for( typename SEQ::data_type::const_iterator i = s.data.begin(); i != s.data.end(); ++i )
      if( !(i->masked()) ) return false;
   return true;
}

inline std::pair<dna_sequence_region_ptr,dna_sequence_region_ptr> make_region_pair( dna_sequence_ptr a, dna_sequence_ptr b, fragment &f ) {
	dna_sequence_region_ptr region_a(new dna_sequence_region(a,f[0]));
	dna_sequence_region_ptr region_b(new dna_sequence_region(b,f[1]));
	return std::make_pair(region_a, region_b);
}

inline int extract_number_of_models( std::istream &infile ) {
	for( int i = 0; i < 7; ++i ) read<double>(infile);
	return read<int>(infile);
}

template<typename VARMAP, typename STR>
void require_option( VARMAP &vm, STR name, int error_code ) {
	using namespace std;
	if( !vm.count(name) ) {
		cerr << "Required parameter --" << name << " missing." << endl;
		exit(error_code);
	}
}

template<typename VARMAP, typename STR>
void require_one_option( VARMAP &vm, STR name_a, STR name_b, int error_code ) {
	using namespace std;
   if( vm.count(name_a) && vm.count(name_b) ) {
      cerr << "Both --" << name_a << " and --" << name_b << " specified. Please use only one." << endl;
      exit(error_code);
   }
	if( !vm.count(name_a) && !vm.count(name_b) ) {
		cerr << "Either --" << name_a << " or --" << name_b << " is required." << endl;
		exit(error_code);
	}
}

template<typename SEQ>
void require_size_above( SEQ &s, unsigned int n ) {
   using namespace std;
   if( s.data.size() <= n ) {
      cerr << "Length of " << s.tags["accession"] << " must be at least " << (n+1) << ".";
      throw runtime_error("Sequence too short.");
   }
}
#endif
