/*
 *  mdna_iupac_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/mdna_iupac.h"
#include "btl/mdna_iupac_symbol.h"
#include "btl/mdna_symbol.h"

std::vector< std::vector<btl::mdna_iupac_symbol> > btl::mdna_iupac_rsets;

btl::mdna_iupac_symbol::mdna_iupac_symbol() : value(0) {}
btl::mdna_iupac_symbol::mdna_iupac_symbol( char const c ) { value = encode(c); }
btl::mdna_iupac_symbol::mdna_iupac_symbol( btl::mdna_iupac_symbol const &o ) : value(o.value) {}

btl::mdna_iupac_symbol::mdna_iupac_symbol( btl::mdna_symbol const &o )  {
	switch(o.value&127) {
		case 0: value = 7; break;
		case 1: value = 3; break;
		case 2: value = 1; break;
		case 4: value = 0; break;
	}
	value |= o.value&128;
}

void btl::mdna_iupac_construct_rsets() {
	using namespace btl;
	mdna_iupac_rsets.resize(15*2);
	mdna_iupac_rsets[mdna_iupac::G.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::C.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::S.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::S.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::T.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::K.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::K.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::Y.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::Y.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::B.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::B.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::B.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::A.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::R.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::R.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::M.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::M.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::V.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::V.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::V.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::W.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::W.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::D.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::D.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::D.index()].push_back(mdna_iupac::G);
	mdna_iupac_rsets[mdna_iupac::H.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::H.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::H.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::N.index()].push_back(mdna_iupac::A);
	mdna_iupac_rsets[mdna_iupac::N.index()].push_back(mdna_iupac::T);
	mdna_iupac_rsets[mdna_iupac::N.index()].push_back(mdna_iupac::C);
	mdna_iupac_rsets[mdna_iupac::N.index()].push_back(mdna_iupac::G);

	// masked version
	mdna_iupac_rsets[15+mdna_iupac::G.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::C.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::S.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::S.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::T.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::K.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::K.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::Y.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::Y.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::B.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::B.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::B.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::A.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::R.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::R.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::M.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::M.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::V.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::V.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::V.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::W.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::W.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::D.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::D.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::D.index()].push_back(mdna_iupac::G.mask());
	mdna_iupac_rsets[15+mdna_iupac::H.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::H.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::H.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::N.index()].push_back(mdna_iupac::A.mask());
	mdna_iupac_rsets[15+mdna_iupac::N.index()].push_back(mdna_iupac::T.mask());
	mdna_iupac_rsets[15+mdna_iupac::N.index()].push_back(mdna_iupac::C.mask());
	mdna_iupac_rsets[15+mdna_iupac::N.index()].push_back(mdna_iupac::G.mask());
}
