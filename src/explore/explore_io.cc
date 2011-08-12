/*
 *  explore_io.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-09-16.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "explore_io.h"
#include <iostream>
#include <fstream>
#include <btl/fasta_reader.h>


ex_sequence_ptr load_sequence( std::string const &file ) {
	using namespace std;
   // load the sequence
   ifstream infile( file.c_str() );
   if( infile.fail() ) throw runtime_error("load_alignment: could not open alignment file");

   ex_sequence_ptr inseq( new ex_sequence() );
   infile >> btl::fasta_reader( *inseq );
   infile.close();

   return inseq;
}

ex_sequence_ptr_vector_ptr load_alignment( std::string const &file ) {
	using namespace std;
   ex_sequence_ptr_vector_ptr myset( new ex_sequence_ptr_vector() );

   // load the sequence
   ifstream infile( file.c_str() );
   if( infile.fail() ) throw runtime_error("load_alignment: could not open alignment file");

   while(1 ) {
		ex_sequence_ptr inseq( new ex_sequence() );
		infile >> btl::fasta_reader( *inseq );
		if( infile.fail() ) break;

      // add to set
      myset->push_back( inseq );
   }
   infile.close();

	return myset;
}
