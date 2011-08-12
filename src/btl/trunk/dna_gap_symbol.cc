/*
 *  dna_gap_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/dna_gap_symbol.h"
#include "btl/mdna_symbol.h"
#include "btl/dna_symbol.h"


btl::dna_gap_symbol::dna_gap_symbol() : value(0) {}
btl::dna_gap_symbol::dna_gap_symbol( char const c ) { value = encode(c); }
btl::dna_gap_symbol::dna_gap_symbol( btl::dna_gap_symbol const &o ) : value(o.value) {}
btl::dna_gap_symbol::dna_gap_symbol( btl::dna_symbol const &o ) : value(o.value) {}
btl::dna_gap_symbol::dna_gap_symbol( btl::mdna_symbol const &o ) : value(o.value&127) {}

