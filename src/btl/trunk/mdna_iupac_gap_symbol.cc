/*
 *  mdna_iupac_gap_symbol.cc
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */

#include "btl/mdna_iupac_gap.h"
#include "btl/mdna_iupac_gap_symbol.h"
#include "btl/mdna_symbol.h"

std::vector< std::vector<btl::mdna_iupac_gap_symbol> > btl::mdna_iupac_gap_rsets;

btl::mdna_iupac_gap_symbol::mdna_iupac_gap_symbol() : value(0) {}
btl::mdna_iupac_gap_symbol::mdna_iupac_gap_symbol( char const c ) { value = encode(c); }
btl::mdna_iupac_gap_symbol::mdna_iupac_gap_symbol( btl::mdna_iupac_gap_symbol const &o ) : value(o.value) {}

btl::mdna_iupac_gap_symbol::mdna_iupac_gap_symbol( btl::mdna_symbol const &o )  {
	switch(o.value&127) {
		case 0: value = 8; break;
		case 1: value = 4; break;
		case 2: value = 2; break;
		case 4: value = 1; break;
	}
	value |= o.value&128;
}

void btl::mdna_iupac_gap_construct_rsets() {
	using namespace btl;
	mdna_iupac_gap_rsets.resize(16*2);
	mdna_iupac_gap_rsets[mdna_iupac_gap::GAP.index()].push_back(mdna_iupac_gap::GAP);
	mdna_iupac_gap_rsets[mdna_iupac_gap::G.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::C.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::S.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::S.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::T.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::K.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::K.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::Y.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::Y.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::A.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::R.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::R.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::M.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::M.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::W.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::W.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::G);
	mdna_iupac_gap_rsets[mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::A);
	mdna_iupac_gap_rsets[mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::T);
	mdna_iupac_gap_rsets[mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::C);
	mdna_iupac_gap_rsets[mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::G);

	// masked version
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::GAP.index()].push_back(mdna_iupac_gap::GAP.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::G.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::C.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::S.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::S.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::T.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::K.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::K.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::Y.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::Y.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::B.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::A.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::R.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::R.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::M.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::M.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::V.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::W.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::W.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::D.index()].push_back(mdna_iupac_gap::G.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::H.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::A.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::T.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::C.mask());
	mdna_iupac_gap_rsets[16+mdna_iupac_gap::N.index()].push_back(mdna_iupac_gap::G.mask());
}
