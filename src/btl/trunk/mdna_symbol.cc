/*
 *  mdna_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/mdna_iupac.h"
#include "btl/mdna_iupac_symbol.h"
#include "btl/mdna_symbol.h"
#include "btl/dna_gap_symbol.h"

btl::mdna_symbol::mdna_symbol() : value(0) {}
btl::mdna_symbol::mdna_symbol( char const c ) { value = encode(c); }
btl::mdna_symbol::mdna_symbol( btl::mdna_symbol const &o ) : value(o.value) {}

// undefined for o being a gap
btl::mdna_symbol::mdna_symbol( btl::dna_gap_symbol const &o ) : value(o.value) {}

// undefined for none mdna characters
btl::mdna_symbol::mdna_symbol( btl::mdna_iupac_symbol const &o ) {
	switch( o.value&127 ) {
		case 7: value = 0; break;
		case 3: value = 1; break;
		case 1: value = 2; break;
		case 0: value = 3; break;
	}
	value |= o.value&128;
}

