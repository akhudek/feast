#ifndef BTL_DNA_GAP__
#define BTL_DNA_GAP__

#include <iostream>
#include <vector>
#include <stdexcept>
#include <boost/array.hpp>
#include <btl/dna_gap_symbol.h>

namespace btl {
	
	class dna_gap {
	public:
		// useful for template constructions
		static int const SIZE = 5;
		
		// reference our symbol type
		typedef dna_gap_symbol symbol;
		
		// list of alphabet symbols
		static boost::array<dna_gap_symbol,5> const alphabet;
		
		// here we encode constant versions of the symbols
		static dna_gap_symbol const A, T, C, G, GAP;
		
   };
}

#endif
