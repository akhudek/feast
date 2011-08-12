/*
 *  dna_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/dna.h"
#include "btl/dna_symbol.h"
#include "btl/mdna_symbol.h"

btl::dna_symbol::dna_symbol() : value(0) {}
btl::dna_symbol::dna_symbol( char const c ) { value = encode(c); }
btl::dna_symbol::dna_symbol( btl::dna_symbol const &o ) : value(o.value) {}

btl::dna_symbol::dna_symbol( btl::mdna_symbol const &o ) : value(o.value&127) {}

// convert from an index to a symbol
btl::dna_symbol btl::dna_symbol::from_index( uint8_t b) {
	using namespace btl;
	switch( b ) {
		case 0: return dna::A;
		case 1: return dna::T;
		case 3: return dna::C;
		case 4: return dna::G;
	}
	throw std::runtime_error( "btl::dna_symbol: invalid index" );
}
