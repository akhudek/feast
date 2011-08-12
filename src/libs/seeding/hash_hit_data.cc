#include "hash_hit_data.h"

// Constructors.
hash_hit_data::hash_hit_data() : match_length(0), num_seqs(0) {}

hash_hit_data::hash_hit_data( unsigned int ml, unsigned int ns ) : match_length(ml), num_seqs(ns) {}


// simple allocator
hash_hit_data_ptr new_hash_hit_data( unsigned int i, unsigned int j ) {
	hash_hit_data_ptr newptr( new hash_hit_data(i,j) );
	return newptr;
}
