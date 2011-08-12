#ifndef MAIN_IO_H__
#define MAIN_IO_H__
// Routines for program input and output
#include <string>
#include <iostream>
#include "dna_sequence.h"

// load a fasta sequence file, disambiguate it, and return a dna_sequence
dna_sequence_ptr load_sequence( std::string const &file );

// Load a set of sequences. Sequences are disambiguated. Limit is the maximum number of sequences read.
dna_sequence_ptr_vector_ptr load_sequence_set( std::string const &file, unsigned int limit = 0 );

#endif
