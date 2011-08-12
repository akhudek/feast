#include "rda_functions.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <btl/algorithm.h>
#include "main_io.h"
#include "nw_model_parameters.h"

nw_model_parameters_ptr load_parameters_from_file( std::string filename ) {
   using namespace std;
   nw_model_parameters_ptr mp( new nw_model_parameters() );

   // open parameter file
   ifstream infile( filename.c_str() );
   if( infile.fail() ) throw runtime_error("set_parameters_from_file: could not open parameter file");
	infile >> *mp;
   infile.close();

   return mp;
}


// the defaults are similar to 2/-3, -5/-2 (NCBI BLAST defaults)
nw_model_parameters_ptr estimate_parameters( dna_sequence_ptr_vector &sequences, double p_d, double p_a, double p_L ) {
   using namespace std;
   nw_model_parameters_ptr mp( new nw_model_parameters());

   blitz::TinyVector<unsigned int,(int)dna_alpha::SIZE> count((unsigned int)0);
   for( dna_sequence_ptr_vector_citer i = sequences.begin(); i != sequences.end(); i++ ) {
      count += symbol_frequencies( **i );
   }

   unsigned int total = sum(count);
   mp->q = (double)count((int)dna_alpha::A)/(double)total, 
           (double)count((int)dna_alpha::T)/(double)total, 
           (double)count((int)dna_alpha::C)/(double)total, 
           (double)count((int)dna_alpha::G)/(double)total; 
   mp->mean_gap_length = p_L;
   mp->pr_open = p_a;
   
   mp->set_p_hky( p_d, 1.0 );
	return mp;
}
