#include "main_io.h"

#include <cctype>
#include <fstream>
#include <stdexcept>
#include <btl/fasta_reader.h>
#include <btl/algorithm.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <algorithm>

#include "entropy.h"
#include "input_sequence.h"

dna_sequence_ptr load_sequence( std::string const &file ) {
   dna_sequence_ptr_vector_ptr myset = load_sequence_set( file, 1 );
   if( myset->size() != 1 ) {
      throw std::runtime_error( "load_sequence: could not read sequence from file" );
   }
   return myset->front();
}

dna_sequence_ptr_vector_ptr load_sequence_set( std::string const &file, unsigned int limit ) {
   using namespace std;
   dna_sequence_ptr_vector_ptr myset = new_dna_sequence_ptr_vector();

   // load the sequence
   ifstream infile( file.c_str() );
   if( infile.fail() ) throw runtime_error("load_sequence: could not open sequence file");

    // random number stream for IUPAC disambiguation
   typedef boost::mt19937 random_number_generator;
   typedef boost::uniform_01<random_number_generator> uniform_random_number_generator;

   // initialize random number generatate for dev random
   entropy random;
   unsigned int seed = (unsigned int)random.get_int();
   if( seed == 0 ) seed = 1; // bug in old boost doesn't allow 0 as a seed
   random_number_generator::result_type rseed(seed);
   random_number_generator rng( rseed );
   uniform_random_number_generator urng( rng );

   bool enable_limit = limit > 0;

   input_sequence inseq;
   while( !enable_limit || limit > 0 ) {
		infile >> btl::fasta_reader( inseq );
		if( infile.fail() ) break;

      limit--;

      // disambiguate IUPAC symbols randomly storing in dna_sequence
      dna_sequence_ptr newseq( new dna_sequence() );
		transform( inseq.data.begin(), inseq.data.end(), back_inserter(newseq->data),
         btl::disambiguate_iupac< random_number_generator, input_alpha, dna_alpha >(rng) );

      // copy over tags
      newseq->tags = inseq.tags;

      // add to set
      myset->push_back( newseq );
   }
   infile.close();
   return myset;
}
