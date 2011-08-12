/*
 *  dna_iupac_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/dna_iupac.h"
#include "btl/dna_iupac_symbol.h"
#include "btl/mdna_symbol.h"

std::vector< std::vector<btl::dna_iupac_symbol> > btl::dna_iupac_rsets;

btl::dna_iupac_symbol::dna_iupac_symbol() : value(0) {}
btl::dna_iupac_symbol::dna_iupac_symbol( char const c ) { value = encode(c); }
btl::dna_iupac_symbol::dna_iupac_symbol( btl::dna_iupac_symbol const &o ) : value(o.value) {}

btl::dna_iupac_symbol::dna_iupac_symbol( btl::mdna_symbol const &o )  {
	switch(o.value&127) {
		case 0: value = 7; break;
		case 1: value = 3; break;
		case 2: value = 1; break;
		case 4: value = 0; break;
	}
}

void btl::dna_iupac_construct_rsets() {
	using namespace btl;
	dna_iupac_rsets.resize(15);
	dna_iupac_rsets[dna_iupac::G.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::C.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::S.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::S.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::T.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::K.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::K.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::Y.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::Y.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::B.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::B.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::B.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::A.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::R.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::R.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::M.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::M.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::V.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::V.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::V.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::W.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::W.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::D.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::D.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::D.index()].push_back(dna_iupac::G);
	dna_iupac_rsets[dna_iupac::H.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::H.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::H.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::N.index()].push_back(dna_iupac::A);
	dna_iupac_rsets[dna_iupac::N.index()].push_back(dna_iupac::T);
	dna_iupac_rsets[dna_iupac::N.index()].push_back(dna_iupac::C);
	dna_iupac_rsets[dna_iupac::N.index()].push_back(dna_iupac::G);
}
